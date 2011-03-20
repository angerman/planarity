(ns geometry
  (:use [incanter.core :only ($= trans matrix bind-rows solve abs to-list sqrt)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Geometry based on a graph structure.

;;------------------------------------------------------------------------------
;; Protocols and Types
(defprotocol GeometryNotation
  (offset        [geom])
  (set-offset    [geom coords] "Sets the geometry offset to coords")
  (reset-offset  [geom] "Sets the geometry offset to nil.")

  (shrink-counter [geom] "Sets the counter to the largest element in nodes.")
  
  (nodes         [geom])
  (add-node      [geom] "Adds a new node to the geometry and returns the node id.")
  (remove-node   [geom node-id] "Removes the node with the given id.")

  (vertices   [geom])
  (add-vertex [geom coords])
  (set-vertex [geom node-id coords])
  (get-vertex [geom node] "Returns the vertex location for a node (int).")
  (remove-vertex [geom node-id])

  (edges      [geom])
  (add-edge   [geom e])
  (remove-edge [geom [a b]])

  (faces    [geom])
  (add-face [geom face])
  (remove-face [geom face]))

(deftype Geometry [^{:volatile-mutable true} nodes
                   ^{:volatile-mutable true} vertices
                   ^{:volatile-mutable true} edges
                   ^{:volatile-mutable true} faces
                   ^{:volatile-mutable true} offset
                   node-counter]
  GeometryNotation
  (offset     [_] offset)
  (set-offset [_ coords] (set! offset coords) coords)
  (reset-offset [_] (set! offset nil))

  (shrink-counter [_] (swap! node-counter (fn [_] (apply max nodes))))
  
  (nodes      [_] nodes)
  (add-node   [_] (let [node (swap! node-counter inc)]
                    (set! nodes (conj nodes node))
                    node))
  (remove-node [_ node] (set! nodes (disj nodes node)))
  
  (vertices [_] vertices)
  (add-vertex [geom coords] (.set-vertex geom (.add-node geom) coords))
  (set-vertex [_ node coords]
              (set! vertices
                    (assoc vertices node (if offset ($= coords + offset) coords)))
              node)
  (get-vertex [_ node] (if offset
                         ($= (get vertices node) - offset)
                         (get vertices node)))
  (remove-vertex [_ node] (set! vertices (dissoc vertices node)))
  
  (edges    [_] edges)
  (add-edge [_ e] (set! edges (conj edges e)) e)
  (remove-edge [_ [a b]] (set! edges (disj edges `(~a ~b))))
  
  (faces    [_] faces)
  (add-face [_ face] (set! faces (conj faces face)) face)
  (remove-face [_ face] (set! faces (disj faces face))))

;;------------------------------------------------------------------------------
;; basic geometry fns
(defn create-geometry []
  (Geometry. #{} {} #{} #{} nil (atom 0)))


;; this is somewhat evil(!) because it binds
;; allows to bind a backup-geometry if
;; *use-backup-geometry* is true.
(def *use-backup-geometry* false)
(if *use-backup-geometry*
  (def *geometry* (if (bound? #'*geometry*)
                    *geometry*
                    (create-geometry)))
  (def *geometry*))

(defn reset-geom []
  (if *use-backup-geometry*
    (def *geometry* (create-geometry))))

(defmacro with-geometry [geometry & body]
  `(binding [*geometry* ~geometry]
     ~@body))

;;------------------------------------------------------------------------------
;; offset
(defn get-offset [] (.offset *geometry*))

(defn set-offset
  ([v] (.set-offset *geometry* v))
  ([a & rest] (set-offset (cons a rest))))

(defn reset-offset [] (.reset-offset *geometry*))

;;------------------------------------------------------------------------------
;; nodes
(defn get-nodes [] (.nodes *geometry*))
(defn add-node [] (.add-node *geometry*))
(defn remove-node [node] (.remove-node *geometry* node))

;;------------------------------------------------------------------------------
;; vertices
(defn get-vertices [] (.vertices *geometry*))
(defn add-vertex
  ([coords] (.add-vertex *geometry* coords))
  ([a & rest] (add-vertex (cons a rest))))

(defn set-vertex
  ([node coords] (.set-vertex *geometry* node coords))
  ([node a & rest] (set-vertex node (cons a rest))))

(defn get-vertex [node] (.get-vertex *geometry* node))

(defn remove-vertex [node] (.remove-vertex *geometry* node))

;;------------------------------------------------------------------------------
;; edges
(defn get-edges [] (.edges *geometry*))
(defn add-directed-edge
  ([e] (let [[a b] e] (.add-edge *geometry* (with-meta `(~a ~b) (meta e)))))
  ([a b] (add-directed-edge [a b])))

(defn add-edge
  ([e] (let [[a b] e] (add-directed-edge (with-meta `(~(min a b) ~(max a b)) (meta e)))))
  ([a b] (add-edge [a b])))

(defn remove-directed-edge
  ([[a b]] (.remove-edge *geometry* `(~a ~b)))
  ([a b] (remove-directed-edge [a b])))

(defn remove-edge
  ([[a b]] (remove-directed-edge `(~(min a b) ~(max a b))))
  ([a b] (remove-edge [a b])))

;;------------------------------------------------------------------------------
;; faces
(defn get-faces [] (.faces *geometry*))
(defn add-face
  ([face] (.add-face *geometry* face))
  ([a & rest] (add-face (cons a rest))))

(defn remove-face
  ([face] (.remove-face *geometry* face))
  ([a & rest] (remove-face (cons a rest))))

;;------------------------------------------------------------------------------
;; helper fn
(defn generate-edges-from-face [face]
  (doseq [[a b] (partition 2 1 [(first face)] face)]
    (add-edge a b))) ;; actually this would have to be
                     ;; `add-directed-edge` but then there would be
                     ;; too many edges

(defn generate-edges-from-faces []
  (doseq [face (.faces *geometry*)]
    (generate-edges-from-face face)))

;;------------------------------------------------------------------------------
;; Interop. Geometry to double arrays
(defn vertex-array []
  (let [nodes (seq (.nodes *geometry*))
        #^"[[D" arr (make-array Double/TYPE (count nodes) (count (get-vertex (first nodes))))]
    (dotimes [i (count nodes)]
      (aset arr (int i) #^doubles (double-array (get-vertex (nth nodes i)))))
    arr))

(defn edge-array []
  (let [edges (seq (.edges *geometry*))
        nodes      (.nodes *geometry*)
        #^"[[I" arr (make-array Integer/TYPE (count edges) 2)
        node-map (apply hash-map (interleave nodes (range (count nodes))))]
    (dotimes [i (count edges)]
      (aset arr (int i) #^ints (int-array (map #(get node-map %) (nth edges i)))))
    arr))

(defn face-array []
  (let [faces (seq (.faces *geometry*))
        nodes      (.nodes *geometry*)
        #^"[[I" arr (make-array Integer/TYPE (count faces) 0)
        node-map (apply hash-map (interleave nodes (range (count nodes))))]
    (dotimes [i (count faces)]
      (aset arr (int i) #^ints (int-array (map #(get node-map %) (nth faces i)))))
    arr))

;;------------------------------------------------------------------------------
;; jReality interop
(defn into-faceset-factory [fsf]
  (.setVertexCount fsf (count (.vertices *geometry*)))
  (.setVertexCoordinates fsf (vertex-array))
  (.setEdgeCount fsf (count (.edges *geometry*)))
  (.setEdgeIndices fsf (edge-array))
  (.setFaceCount fsf (count (.faces *geometry*)))
  (.setFaceIndices fsf (face-array))
  (.update fsf))

(defn into-edgeset-factory [esf]
  (.setVertexCount esf (count (.vertices *geometry*)))
  (.setVertexCoordinates esf (vertex-array))
  (.setLineCount esf (count (.edges *geometry*)))
  (.setEdgeIndices esf (edge-array))
  (.update esf))

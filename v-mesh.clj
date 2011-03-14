;; -*- coding: utf-8 -*-
(ns v-mesh
  (:import [java.awt Color])
  (:import [java.awt.event KeyListener KeyEvent])
  (:import [de.jreality.shader.ShaderUtility])
  (:import [[de.jreality.scene Appearance Geometry SceneGraphComponent]])
  (:import [[de.jreality.scene.data Attribute StorageModel IntArrayArray$Array]])
  (:import [[de.jreality.tools DragEventTool]])
  (:import [[de.jreality.plugin JRViewer]])
  (:import [[de.jreality.geometry IndexedFaceSetFactory QuadMeshFactory]])
  (:import [de.jreality.geometry.QuadMeshFactory])
  (:use [incanter.core])
  (:use [clojure.contrib.seq-utils :only (indexed)])) ;; For $=, decompose-eig and general matrix
;; computation

(def id [[1 0 0 0]
         [0 1 0 0]
         [0 0 1 0]
         [0 0 0 1]])

(def origin [1 0 0 0])

(defn rot-x [t]
      (matrix [[1 0 0 0] [0 1 0 0] [0 0 (cos t) (* -1 (sin t))] [0 0 (sin t) (cos t)]]))
(defn rot-z [a]
      (matrix [[1 0 0 0] [0 (cos a) (* -1 (sin a)) 0] [0 (sin a) (cos a) 0] [0 0 0 1]]))
(defn translate [[x y z]]
      (matrix [[1 0 0 0] [x 1 0 0] [y 0 1 0] [z 0 0 1]]))
(defn move [t a l]
  ($= (rot-x t) <*> (rot-z a) <*> (translate [l 0 0])))

(defn norm [[x y z]]
  (sqrt ($= x ** 2 + y ** 2 + z ** 2)))

(defn rot [[x y z] [a b c] [u v w] t]
  (let [ct (cos t)
        st (sin t)
        ux ($= u * x + v * y + w * z)
        denom ($= u ** 2 + v ** 2 + w ** 2)
        n (sqrt denom)]
    [($= (a * (v ** 2 + w ** 2) + u * (ux - b * v - c * w) + ((x - a) * (v ** 2 + w ** 2) + u * (b * v + c * w - v * y - w * z)) * ct + n * (b * w - c * v - w * y + v * z) * st) / denom)
     ($= (v * (u ** 2 + w ** 2) + v * (ux - a * u - c * w) + ((y - b) * (u ** 2 + w ** 2) + v * (a * u + c * w - u * x - w * z)) * ct + n * (c * u - a * w + w * x - u * z) * st) / denom)
     ($= (c * (u ** 2 + v ** 2) + w + (ux - a * u - b * v) + ((z - c) * (u ** 2 + v ** 2) + w * (a * u + b * v - u * x - v * y)) * ct + n * (a * v - b * u - v * x + u * y) * st) / denom)
     1]))

(defn dehomogenize [v]
  (let [v (vec (reverse v))]
    (vec (reverse ($= (rest v) / (first v))))))

(def pi:8 (/ Math/PI 32))

(def u-spec '((0.2 0.3 1)
              (0.3 0.4 1)
              (0.4 0.5 1)
              (0.5 0.4 1)
              (0.4 0.3 1)
              (0.3 0.2 1)))
(def v-spec '((0.2 0.3 1)
              (0.3 0.4 1)
              (0.4 0.5 1)
              (0.5 0.4 1)
              (0.4 0.3 1)
              (0.3 0.2 1)))
;;(def u-spec (take 12 (repeat (list 0 0 1))))
;;(def v-spec (take 12 (repeat (list 0 0 1))))



(def u-points (points-from-spec u-spec))
(def v-points (map #($= (rot-z (* 2 (/ Math/PI 4))) <*> %) (points-from-spec v-spec)))

(def u-points (map #($= (translate [0 0 2]) <*> %) u-points))
(def v-points (map #($= (translate [0 0 2]) <*> %) v-points))
(def origin ($= (translate [0 0 2]) <*> origin))


;; Plücker Coordinates
(defn pdet [x y i j]
  "Computes the determinant of rows, i,j of [x y]"
  (- (* (nth x i) (nth y j)) (* (nth x j) (nth y i))))

(defn pcross [x y]
  [(pdet x y 1 2)
   (pdet x y 2 0)
   (pdet x y 0 1)])

(defn pline [x y]
  "Assumes x and y in homogenous coordinates. Result is
   [p_01,p_02,p_03,p_23,p_21,_p12"
  [(pdet x y 0 1)
   (pdet x y 0 2)
   (pdet x y 0 3)
   (pdet x y 2 3)
   (pdet x y 3 1)
   (pdet x y 1 2)])

(defn pline-meet [[d0 d1 d2 m0 m1 m2] [e0 e1 e2 n0 n1 n2]]
  (vec (cons (+ (* d0 n0) (* d1 n1) (* d2 n2)) (pcross [m0 m1 m2] [n0 n1 n2]))))

(defn line-meet [x dx y dy]
  (let [a (pcross dx dy)
        b (pcross ($= y - x) dy)
        scale (if-not (= (nth b 0) 0)
                (/ (nth a 0) (nth b 0))
                (if-not (= (nth b 1) 0)
                  (/ (nth a 1) (nth b 1))
                  (if-not (= (nth b 2) 0)
                    (/ (nth a 2) (nth b 2)))))]
    ($= x + scale * dx)))

;; Working with a mesh
(def mesh
  (partition 7 (take (* 7 7) (repeatedly #(ref nil)))))

(defn node [u v]
  (nth (nth mesh u nil) v nil))

(defn set-node-pos [u v pos]
  (dosync
   (if-let [n (node u v)]
     (ref-set n pos))))

(declare node-du node-dv)

(defn node-pos [u v]
  (if-let [n (node u v)]
    (if @n
      @n
      (do ;; compute node
        (let [pu (node-pos (dec u) v)
              dv (node-dv  (dec u) v)
              pv (node-pos u (dec v))
              du (node-du  u (dec v))]
          (if (and pu du pv dv)
            (let [pos (cons 1 (line-meet (rest pu) (rest du) (rest pv) (rest dv)))]
              (set-node-pos u v pos)
              pos)))))))

(defn node-du [u v]
  (let [n  (node-pos u v)
        pn (node-pos (dec u) v)]
    (if (and n pn)
      ($= n - pn))))

(defn node+du [u v]
  (let [n  (node-pos u v)
        nn (node-pos (inc u) v)]
    (if (and n nn)
      ($= n - nn))))

(defn node-dv [u v]
  (let [n (node-pos u v)
        pn (node-pos u (dec v))]
    (if (and n pn)
      ($= n - pn))))

(defn node+dv [u v]
  (let [n  (node-pos u v)
        nn (node-pos u (inc v))]
    (if (and n nn)
      ($= n - nn))))

(defn deho [v]
  ($= (rest v) / (first v)))

(defn node-pos-double-array [u v]
  (if-let [n (node-pos u v)]
    (double-array ($= (rest n) / (first n)))))

;; Setup Mesh
(defn setup-mesh []
  (set-node-pos 0 0 origin)
  (doseq [i (range (length u-points))]
    (set-node-pos (inc i) 0 (nth u-points i)))
  (doseq [i (range (length v-points))]
    (set-node-pos 0 (inc i) (nth v-points i))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; local planar deform construction
(defn neg [x] (* -1 x))
(defn skew-mat [[x y z]]
  (matrix [[ 0 (neg z)  y]
           [ z  0 (neg x)]
           [(neg y)  x  0]]))

(declare three-to-one)
(defn dual-node-pos [u v]
  (three-to-one [(node-pos u (dec v))
                 (node-pos u v)
                 (node-pos (dec u) v)]))

(defn dual-node-pos-double-array [u v]
  (if-let [n (dual-node-pos u v)]
    (double-array ($= (rest n) / (first n)))))

;;      
;;     L
;;   G H I
;; M D E F K
;;   A B C
;;     J
;;
(defn e3 [n]
  (vec
   (for [i (range 3)]
     (if (= i n) 1 0))))
(defn dec2 [n]
  (dec (dec n)))
(defn inc2 [n]
  (inc (inc n)))

(declare psf qmf2)

(defn update-point [index pos]
  (let [ps (. psf getPointSet)
        n  (. ps getNumPoints)
        pts (make-array Double/TYPE n 3)
        coords de.jreality.scene.data.Attribute/COORDINATES
        sm (. de.jreality.scene.data.StorageModel/DOUBLE_ARRAY array 3)]
    (. (. ps getVertexAttributes coords) toDoubleArrayArray pts)
    (aset pts index pos)
    (. ps setVertexAttributes coords (. sm createReadOnly pts))))

(defn render-plane [E q1 Q1 Q2 Q3 Q4]
  (let [;; first compute the four points
        q2 ($= Q1 <*> q1)
        q3 ($= Q2 <*> q2)
        q4 ($= Q3 <*> q3)]
        ;; normalize them to lie in the plane of E
    (map #(cons 1.0 ($= -1 * ( % / ( (trans (rest E)) <*> %) ) ))
         (list q1 q2 q3 q4))))

(defn set-pts [pts offset]
  (doseq [[idx pt] (indexed pts)]
    (update-point (+ offset idx)
                  (double-array (rest pt)))))

(defn three-to-one [pts]
  (let [lhs (apply bind-rows (map (comp trans matrix rest) pts))
        rhs (matrix (map (comp #(* -1 %) first) pts))]
    (cons 1.0 (solve lhs rhs))))

(defn vec-length [v]
  (sqrt ($= (trans v) <*> v)))

(defn normalize [v]
  ($= v / (vec-length v)))

(defn compute-points [u v q1 q2 q3 q4]
  (let [b (three-to-one [(node-pos (dec u) (dec v))
                         (node-pos u (dec v))
                         (node-pos u (dec2 v))])
        f (three-to-one [(node-pos (inc u) (dec v))
                         (node-pos (inc u) v)
                         (node-pos (inc2 u) v)])
        h (three-to-one [(node-pos (inc u) (inc v))
                         (node-pos u (inc v))
                         (node-pos u (inc2 v))])
        d (three-to-one [(node-pos (dec u) (dec v))
                         (node-pos (dec u) v)
                         (node-pos (dec2 u) v)])]
    ;; compute new points
    (list (three-to-one [b q1 q2]) ;; new b
          (three-to-one [f q2 q3]) ;; new f
          (three-to-one [h q3 q4]) ;; new h
          (three-to-one [d q4 q1])))) ;; new d

(defn deform [u v p]
  "computes the deformation positions for point (u,v) with new coordinates of p"
  (let [Me (skew-mat (normalize (take 3 p)))
        Mj (skew-mat (normalize (rest (node-pos u (dec2 v)))))
        Mb (skew-mat (normalize (rest (node-pos u (dec  v)))))
        Mc (skew-mat (normalize (rest (node-pos (dec u) (dec v)))))
        Mk (skew-mat (normalize (rest (node-pos (inc2 u) v))))
        Mf (skew-mat (normalize (rest (node-pos (inc u) v))))
        Mi (skew-mat (normalize (rest (node-pos (inc u) (inc v)))))
        Ml (skew-mat (normalize (rest (node-pos u (inc2 v)))))
        Mh (skew-mat (normalize (rest (node-pos u (inc v)))))
        Mg (skew-mat (normalize (rest (node-pos (inc u) (dec v)))))
        Md (skew-mat (normalize (rest (node-pos (dec u) v))))
        Mm (skew-mat (normalize (rest (node-pos (dec2 u) v))))
        Ma (skew-mat (normalize (rest (node-pos (dec u) (dec v)))))
        ;; Q1: q1 -> q2
        [Q4 Q3 Q2 Q1] (map (fn [[a b c]] ($= a <*> b <*> c))
                           (partition 3 
                                      (list Ma Md Mm
                                            Mg Mh Ml
                                            Mi Mf Mk
                                            Mc Mb Mj)))
        M (reduce #($= % <*> %2) (list Q4 Q3 Q2 Q1))
        ;; find the eigenvalues
        ;; and select the eigenvector corresponding
        ;; to the larger eigenvalue
        p0 (cons 1.0 (take 3 p))
        eig (decomp-eigenvalue M)
        val (map abs (:values eig))
        _ (prn (:values eig))
        _ (prn (:vectors eig))
        min-val (apply min val)
        evs (for [[k v] (indexed val) :when (not= v min-val)]        
              k)]
    (doseq [[i k] (indexed evs)]
      (let [[q1 q2 q3 q4] (render-plane p0 ($= (:vectors eig) <*> (e3 k)) Q1 Q2 Q3 Q4)]
        (set-pts (compute-points u v q1 q2 q3 q4) (inc (* 4 i)))))
    (update-point 0 (double-array (take 3 p)))
    (let [pt (reduce #($= % + %2)
                     (for [du [-1 0 1]
                           dv [-1 0 1] :when (not= du dv 0)]
                       (rest (node-pos (+ u du) (+ v dv)))))]
      (update-point 9 (double-array ($= pt / 8)))
      (let [pts (reduce concat
                    (map (comp to-list)
                         (for [du (list -1 0 1 2)
                               dv (list -1 0 1 2)]
                           ($= (rest (dual-node-pos (+ u du) (+ v dv))) - pt))))]
        (doto qmf2
          (.setVertexCoordinates
           (double-array pts))
          (.update))))
    
;;    (set-pts [b f h d] 9)
    ))



;; jReality Interop
(defn node-array []
  (let [#^"[[[D" arr (make-array Double/TYPE (count mesh) (count (first mesh)) 3)]
    (dotimes [u (count mesh)]
      (dotimes [v (count (nth mesh u))]
        (aset arr (int u) (int v) #^doubles (node-pos-double-array u v))))
    arr))

(defn dual-node-array []
  (let [#^"[[[D" arr (make-array Double/TYPE (dec (count mesh)) (dec (count (first mesh))) 3)]
    (dotimes [u (dec (count mesh))]
      (dotimes [v (dec (count (nth mesh u)))]
        (aset arr (int u) (int v) #^doubles (dual-node-pos-double-array u v))))
    arr))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; QuadMeshFactory
(setup-mesh)

(prn (three-to-one [[1 0 0 1] [1 0 1 1] [1 1 1 1]]))

#_(prn (bind-columns (node-pos (dec 2) 1)
                   (node-dv (dec 2) 1)
                   (node-pos 2 (dec 1))
                   (node-du 2 (dec 1))))


#_(prn (apply bind-columns (for [u (list 0 1 2 3 4 5)] (matrix (dual-node-pos u 0)))))

;;(prn (count mesh))
;; (dotimes [u (count mesh)]
;;   (dotimes [v (count mesh)]
;;     (prn [(aget (node-array) u v 0)
;;           (aget (node-array) u v 1)
;;           (aget (node-array) u v 2)])))

(def qmf (de.jreality.geometry.QuadMeshFactory.))
(doto qmf
  (.setVLineCount (count mesh))
  (.setULineCount (count mesh))
  (.setClosedInUDirection false)
  (.setClosedInVDirection false)
  (.setVertexCoordinates (node-array))
  (.setGenerateFaceNormals true)
  (.setGenerateTextureCoordinates true)
  (.setGenerateEdgesFromFaces true)
  (.setEdgeFromQuadMesh true)
  (.update))

(def qmf2 (de.jreality.geometry.QuadMeshFactory.))
(doto qmf2
  (.setVLineCount 4)
  (.setULineCount 4)
  (.setClosedInUDirection false)
  (.setClosedInVDirection false)
  (.setVertexCoordinates (double-array (take (* 4 4 3) (repeat 0))))
  (.setGenerateFaceNormals true)
  (.setGenerateTextureCoordinates true)
  (.setGenerateEdgesFromFaces true)
  (.setEdgeFromQuadMesh true)
  (.update))


(def psf (de.jreality.geometry.IndexedLineSetFactory.))

(doto psf
  (.setVertexCount 13)
  (.setVertexCoordinates (double-array (take (* 13 3) (repeat 0))))
  (.setEdgeCount 12)
  (.setEdgeIndices (int-array [0 1 0 2 0 3 0 4
                               0 5 0 6 0 7 0 8
                               1 5 2 6 3 7 4 8]) 2)
  (.setEdgeColors (let [arr (make-array java.awt.Color 12)]
                    (doseq [i (range 12)]                     
                      (aset arr i (if (< i 4) java.awt.Color/green
                                      (if (< i 8) java.awt.Color/orange
                                          java.awt.Color/black))))
                    arr))
;;  (.setGenerateFaceNormals true)
  (.update))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; DragTool
(def dt (de.jreality.tools.DragEventTool.))
(defn movable [u v]
  (and (> u 1) (> v 1)
       (< u 5) (< v 5)))

(.addPointDragListener dt (proxy [de.jreality.tools.PointDragListener] []
                            (pointDragStart [e])
                            (pointDragEnd [e])
                            (pointDragged [e]
                                          (let [idx (.getIndex e)
                                                u (mod idx 7)
                                                v (/ (- idx u) 7)]
;;                                            (prn (list u v))
                                            (if (movable u v)
                                              (let [old (node-pos u v)
                                                    old (rest old)
                                                    new (map #(aget (.getPosition e) %)
                                                             (range (alength (.getPosition e))))
                                                    new (take 3 new)]
                                                (deform v u new)))))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Appearance setup
(defn setup-appearance [ap]
  (let [gs (de.jreality.shader.ShaderUtility/createDefaultGeometryShader ap true)
        ps (. gs createPolygonShader "default")
        ls (. gs createLineShader "default")
        pts (. gs createPointShader "default")]
    (doto gs
      (.setShowFaces true)
      (.setShowLines true)
      (.setShowPoints true))
    (doto ps
      (.setDiffuseColor Color/blue))
    (doto ls
      (.setDiffuseColor Color/yellow)
      (.setTubeRadius 0.03))
    (doto pts
      (.setDiffuseColor Color/red)
      (.setPointRadius 0.05)))
  ap)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Component Graph assembly and Viewer initialization
(def cmp (de.jreality.scene.SceneGraphComponent.))
(def dual-cmp (de.jreality.scene.SceneGraphComponent.))
(def subcmp (de.jreality.scene.SceneGraphComponent.))
(doto subcmp
  (.setGeometry (.getGeometry psf)))
(doto dual-cmp
  (.setGeometry (.getGeometry qmf2))
  (.setAppearance (setup-appearance (de.jreality.scene.Appearance.))))
(doto cmp
  (.setGeometry (.getGeometry qmf))
  (.setAppearance (setup-appearance (de.jreality.scene.Appearance.)))
  (.addTool dt)
  (.addChild subcmp)
  (.addChild dual-cmp))


(def viewer (de.jreality.plugin.JRViewer/display cmp))
;;(def dual-viewer (de.jreality.plugin.JRViewer/display dual-cmp))


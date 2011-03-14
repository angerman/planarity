(ns geometry.subdivision.doo-sabin
  (:use [incanter.core :only ($= bind-columns)])
  (:use [geometry.utils :only (faces-for-edge edge-in-face face-to-directed-edges faces-for-vertex)]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Doo Sabin
(defn doo-sabin-subd-weights [n]
  (cons (+ 1/4 (/ 5 (* 4 n)))
        (for [i (range 1 n)]
          (/ (+ 3 (* 2 (Math/cos (/ (* 2 Math/PI i) n))))
             (* 4 n)))))

(defn doo-sabin-subd-face-points []
  {:fv->v
   (apply merge
          (apply concat
                 (for [face (geometry/get-faces)]
                   (let [n (count face)
                         w (doo-sabin-subd-weights n)]
                     (for [i (range n)]
                       (let [vertices (take n (drop i (cycle face)))
                             vertex   (first vertices)
                             mat      (apply bind-columns (map geometry/get-vertex vertices))]
                         {(list face vertex) (apply geometry/add-vertex ($= mat <*> w))}))))))})

(defn doo-sabin-subd-face-faces [struct]
  (doseq [face (set (map first (keys (:fv->v struct))))]
    (apply geometry/add-face
           (map #(get (:fv->v struct) (list face %)) face)))
  struct)

(defn doo-sabin-subd-edge-faces [struct]
  (let [faces (set (map first (keys (:fv->v struct))))
        fv->v  #(get (:fv->v struct) (list % %2))
        edges (geometry.utils/edges-from-faces faces)]
    (doseq [[a b] edges]
      (let [face-l (or (first (faces-for-edge faces [a b]))
                       (prn "[Edge -> Face] WARNING: No Left face!")
                       (second (faces-for-edge faces [b a])))
            face-r (or (first (faces-for-edge faces [b a]))
                       (prn "[Edge -> Face] WARNING: No Right face!")
                       (second (faces-for-edge faces [a b])))]
        (geometry/add-face
         (fv->v face-l a)
         (fv->v face-l b)
         (fv->v face-r b)
         (fv->v face-r a))))
    struct))

(defn doo-sabin-order-faces [faces vertex face]
  (if face
    (if (empty? faces)
      (list face)
      (let [face-for-edge (fn [e] (first (filter (partial edge-in-face e) faces)))
            edge (first (filter #(= vertex (first %)) (face-to-directed-edges face)))
            edge2 (first (filter #(= vertex (second %)) (face-to-directed-edges face)))
            next-face (or (face-for-edge (reverse edge))
                          (prn "[Vertex -> Face] WARNING: orientation issue!")
                          (face-for-edge edge)
                          (face-for-edge edge2)
                          (face-for-edge (reverse edge2)))]
        (if next-face
          (cons face (doo-sabin-order-faces (disj faces next-face) vertex next-face)))))
    (doo-sabin-order-faces (set (rest faces)) vertex (first faces))))

(defn doo-sabin-subd-vertex-faces [struct]
  (let [faces (set (map first (keys (:fv->v struct))))
        face-for-edge (fn [e] (first (filter (partial edge-in-face e) faces)))
        fv->v #(get (:fv->v struct) (list % %2))
        edges (geometry.utils/edges-from-faces faces)
        vertices (set (apply concat edges))]
    (doseq [vertex vertices]
      (let [ofaces (doo-sabin-order-faces (faces-for-vertex faces vertex) vertex nil)]
        (apply geometry/add-face
               (for [face ofaces]
                 (fv->v face vertex)))))
    struct))

(defn doo-sabin-subd-cleanup [struct]
  (let [faces (set (map first (keys (:fv->v struct))))]
    (doseq [face faces]
      (apply geometry/remove-face face))
    (doseq [[a b] (geometry.utils/edges-from-faces faces)]
      (geometry/remove-edge a b))
    (doseq [vertex (set (apply concat faces))]
      (geometry/remove-vertex vertex)
      (geometry/remove-node vertex))))

(defn doo-sabin-subd []
  (-> (doo-sabin-subd-face-points)
      (doo-sabin-subd-face-faces)
      (doo-sabin-subd-edge-faces)
      (doo-sabin-subd-vertex-faces)
      (doo-sabin-subd-cleanup))
  (geometry/generate-edges-from-faces))

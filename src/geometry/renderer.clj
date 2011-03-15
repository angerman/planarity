(ns geometry.renderer
  (:use [incanter.core :only ($= matrix bind-columns trans)])
  (:use [geometry :only (with-geometry get-nodes get-edges get-vertices get-vertex get-faces)])
  (:use [geometry.utils :only (unit-vec skew-mat norm centroid)])
  (:use [clojure.string :only (join)]))

;;------------------------------------------------------------------------------
;; slightly more advanced renderer (than the previous one)
;;
;; TODO
;; - faces
;; - lighting 

(def *node-counter*)
(defn next-id []
  (swap! *node-counter* inc))

(defmacro with-counter [ & body]
  `(binding [*node-counter* (atom 0)]
     ~@body))

(defrecord Node  [id x y])
(defrecord Point [node depth])
(defrecord Edge  [node-a node-b depth-a depth-b])
(defrecord Face  [nodes depths])

(defprotocol TikZRenderer
  (render [_] "TikZ representation.")
  (depth  [_] "Computes the average depth."))

(extend-protocol TikZRenderer
  Node
  (render [n] (format "\\coordinate (node-%d) at (%+2.4f,%+2.4f);\n" (:id n) (float (:x n)) (float (:y n))))
  Point
  (depth  [p] (:depth p))
  (render [p] (format "\\fill[black!%d] (node-%d) circle (1.5pt);\n" (max 0 (min 100 (Math/round (* 100 (:depth p))))) (:node p)))
  Edge
  (depth  [e] (geometry.utils/avg (:depth-a e) (:depth-b e)))
  (render [e] (format "\\draw[black!%d] (node-%d) -- (node-%d);\n"
                      (max 0 (min 100 (Math/round (* 100 (depth e)))))
                      (:node-a e) (:node-b e)))
  Face
  (depth  [f] (geometry.utils/avg (:depths f)))
  (render [f] (format "\\fill[black!%d] %s -- cycle;\n"
                      (max 0 (min 80 (Math/round (* 80 (depth f)))))
                      (join " -- " (map #(format "(node-%d)" %) (:nodes f))))))

(defn make-node [x y]
  (Node. (next-id) x y))

(defn make-edge [n-a n-b d-a d-b]
  (Edge. n-a n-b d-a d-b))


(defn cross-prod [[a b c] [d e f]]
  ;; this could also be accomplished
  ;; with the skew-mat fn.
  [(- (* b f) (* c e))
   (+ (* a f) (* c d))
   (- (* a e) (* b d))])

(defn rot-mat [x y z]
  (bind-columns [1 0 0 0]
                (cons 0 x)
                (cons 0 y)
                (cons 0 z)))

(defn trans-mat [[x y z]]
  (matrix [[1 0 0 0]
           [x 1 0 0]
           [y 0 1 0]
           [z 0 0 1]]))

(defn scale-mat [[x y z]]
  (matrix [[1 0 0 0]
           [0 x 0 0]
           [0 0 y 0]
           [0 0 0 z]]))

(defn view-matrix [target eye up]
  "Computes the view matrix. Such that the coordinate system lies
   at the eye and points (-z) towards the target. With the y coordinate
   approximating the `up` direction.

   In case of up || (eye - target) behaviour is undefined."
  (let [up    (unit-vec up)
        look  (unit-vec ($= eye - target))
        right (unit-vec (cross-prod up look))
        up    (unit-vec (cross-prod look right))
        R     (trans (rot-mat right up look))
        T     (trans-mat ($= -1 * eye))]
    ($= R <*> T )))

(defn render-geometry [geom eye up near far & {:keys [magnify] :or {magnify 1.0}}]
  "No clipping is performed"
  (with-geometry geom
    (let [target (centroid (vals (get-vertices)))
          view-m ($= (scale-mat [magnify magnify 1.0]) <*> (view-matrix target eye up))
          n->xy  (fn [n] (let [v* ($= view-m <*> (cons 1 (get-vertex n)))]
                          ($= (take 2 (rest v*)) / (-1 * (last v*)))))
          n->depth (fn [n] (let [v*   ($= view-m <*> (cons 1 (get-vertex n)))
                                z    ($= -1 * (last v*))
                                dist (- far near)]
                            (- 1 (/ z dist))))]
      (with-counter
        (let [n->tikz-n (reduce merge (for [n (get-nodes)] {n (apply make-node (n->xy n))}))
              points    (for [n (get-nodes)] (Point. (:id (n->tikz-n n)) (n->depth n)))
              edges     (for [[a b] (get-edges)] (Edge. (:id (n->tikz-n a))
                                                        (:id (n->tikz-n b))
                                                        (n->depth a)
                                                        (n->depth b)))
              faces     (for [face (get-faces)] (Face. (map #(:id (n->tikz-n %)) face) (map n->depth face)))]
          (doall
           (map render (concat (vals n->tikz-n) (sort-by depth (concat points edges faces))))))))))

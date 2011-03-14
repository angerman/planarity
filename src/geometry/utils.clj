(ns geometry.utils
  (:use [incanter.core :only ($= trans sqrt matrix abs solve bind-rows to-list)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)]))

;;------------------------------------------------------------------------------
;; edge and face fn
(defn edge
  ([[a b]] `(~(min a b) ~(max a b)))
  ([a b] (edge [a b])))

(defn face-to-directed-edges [face]
  (partition 2 1 [(first face)] face))

(defn face-to-edges [face]
  (map #(apply edge %)
       (face-to-directed-edges face)))

(defn edges-from-faces [faces]
  (->> (for [face faces]
         (partition 2 1 [(first face)] face))
       (apply concat)
       (map edge)
       (set)))

(defn edge-in-face [edge face]
  (some #{edge} (face-to-directed-edges face)))

(defn edges-for-vertex [edges vertex]
  (set (filter #(some #{vertex} %) edges)))

(defn faces-for-edge [faces edge]
  (filter #(edge-in-face edge %) faces))

(defn faces-for-vertex [faces vertex]
  (set (filter #(some #{vertex} %) faces)))

;;------------------------------------------------------------------------------
;; general tools
(defn avg
  ([lst] (/ (reduce + lst) (count lst)))
  ([a & rest] (avg (cons a rest))))

(defn centroid
  ([pts] (let [N (count pts)] (map #(/ % N) (reduce #(map + % %2) pts))))
  ([p & rest] (centroid (cons p rest))))

(defn norm [v]
  (sqrt ($= (trans v) <*> v)))

(defn unit-vec [v]
  ($= v / (norm v)))

(defn neg [x] (* -1 x))

(defn skew-mat [[x y z]]
  (matrix [[ 0 (neg z)  y]
           [ z  0 (neg x)]
           [(neg y)  x  0]]))

(defn parrallel? [v1 v2]
  (let [uv1 (unit-vec v1)
        uv2 (unit-vec v2)]
    (approx= (abs ($= (trans uv1) <*> uv2)) 1.0 1E-16)))

(defn plane [a b c]
  (try
    (solve (apply bind-rows (map to-list [a b c]))
           (matrix [-1 -1 -1]))
    (catch IllegalArgumentException e
      (unit-vec ($= (skew-mat ($= b - a)) <*> ( b - c))))))

(defn point [p q r] ;; dual version of plane
  (plane p q r))

(defn polygon-length [pts]
  (reduce +
          (for [[a b] (partition 2 1 [(first pts)] pts)]
            (norm ($= b - a)))))

;; the idea is
;; compute the directions of the edges successively
;; and compute their scalar product. if the sign
;; of the scalar product does not change, it's convex.
;; for 3d we will have to fix a plane first.
(defn is-convex? [pts]
  (let [edges (map (fn [[a b]] ($= b - a)) (partition 2 1 [(first pts)] pts))
        plane ($= (skew-mat (first edges)) <*> (second edges))
        cons-edges (partition 2 1 [(first edges)] edges)]
    (reduce #(and % %2)
            (map (fn [[a b]] (pos? ($= (trans plane) <*> (skew-mat a) <*> b))) cons-edges))))

(defn quad-area [[a b c d]]
  (* 1/2
     (norm ($= ( (skew-mat ($= c - a)) <*> ( d - b ))))))

;;------------------------------------------------------------------------------
;; line and plane intersections
(defn line-through-points [a b]
  (list a ($= b - a)))

(defn line-meet [x dx y dy]
  (let [a ($= (skew-mat dx) <*> dy)
        b ($= (skew-mat ($= y - x)) <*> dy)
        scale (if-not (= (nth b 0) 0)
                (/ (nth a 0) (nth b 0))
                (if-not (= (nth b 1) 0)
                  (/ (nth a 1) (nth b 1))
                  (if-not (= (nth b 2) 0)
                    (/ (nth a 2) (nth b 2)))))]
    ($= x + scale * dx)))

(defn line-plane-intersection [[off dir] p]
  (if (< (abs ($= (trans dir) <*> p)) 1E-17)
    (println "Cannot compute line-plane intersection. Line is parrallel to plane!")
    (let [t ($= ( -1 - (trans off) <*> p) / ( (trans dir) <*> p) )]
      ($= off + t * dir))))

(defn plane-plane-intersection [[a b c] [d e f]]
  (try
    (let [[g h i] ($= (skew-mat [a b c]) <*> [d e f]) ;; direction
          off (solve (matrix [[a b c]
                              [d e f]
                              [g h i]])
                     (matrix [-1 -1 0]))]
      ($= (skew-mat off) <*> ( off + [g h i] )))
    (catch java.lang.IllegalArgumentException ex
      (println (str "Cannot compute intersection of given planes with normal vectors: " [a b c] " and " [d e f] " (" (.getMessage ex) ")")))))

(defn parallel-plane-through-point [h p]
  (let [normal (unit-vec h)
        scale  ($= -1 * (trans normal) <*> p)]
    ($= normal / scale)))

(defn point-in-plane [p h]
  "Orhtogonal projectiong of point p in plane h"
  (let [scale ($= ( (trans h) <*> p + 1 ) / ( (trans h) <*> h ) )]
    ($= p - scale * h)))

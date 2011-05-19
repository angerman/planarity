(ns geometry.subdivision.planar
  (:require [tikz :as t])
  (:use [geometry :only (add-vertex get-vertex add-edge remove-edge get-edges)])
  (:use [geometry.utils :only (centroid edges-for-vertex face-to-directed-edges edge quad-area polygon-length
                                        is-convex? parrallel? face-to-edges unit-vec skew-mat norm
                                        plane line-plane-intersection line-through-points planar? planar-dist cross-prod)])
  (:use [geometry.renderer :only (render-geometry)])
  (:use [planar-dual :only (M)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)])
  (:use [incanter.core :only ($= trans bind-columns)])
  (:use [clojure.contrib.duck-streams :only (with-out-writer)])
  (:use [planar.dual-solver :only (solve with-face)])
  (:use [clojure.string :only (join)]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Planar SubD
;; for each edge, add 1/3, 2/3 points
;; adjust old vertices.

(def *ps)
(def *faces)
(def *counter* (atom 0))
(def *use-backup* true)

;; the weights for the vertex repositioning.
(def *vertex-weight* 1/3)
(def *centroid-weight* 2/3)
(def INSERT-A 2/3)
(def INSERT-B 1/3)

(defn planar-subd-edge-points []
  {:vertices (set (reduce concat (geometry/get-edges)))
   :ev->v (reduce merge
                  (for [[a b] (geometry/get-edges)]
                    (let [A (geometry/get-vertex a)
                          B (geometry/get-vertex b)
                          a* (apply geometry/add-vertex ($= INSERT-A * A + INSERT-B * B))
                          b* (apply geometry/add-vertex ($= INSERT-B * A + INSERT-A * B))]
                      {`((~a ~b) ~a) a*
                       `((~a ~b) ~b) b*})))})

(defn planar-subd-adjust-vertices [struct]
  (doseq [vertex (:vertices struct)]
    (let [edges (edges-for-vertex (geometry/get-edges) vertex)]
      (let [vlocation (geometry/get-vertex vertex)
            ecentroid (apply centroid 
                             (map geometry/get-vertex
                                  (map #(get (:ev->v struct) `(~% ~vertex)) edges)))
            n (count edges)]
        (geometry/set-vertex vertex ($=   ( *vertex-weight* * vlocation )
                                          + ( *centroid-weight* * ecentroid )))))))

(defn planar-subd-compute-plane-point-for-face [struct face]
  (let [edges (face-to-directed-edges face)]
    (let [edge-points 
          (for [[a b] edges]
            (list
             (get (:ev->v struct) `(~(edge a b) ~a))
             (get (:ev->v struct) `(~(edge a b) ~b))))
          ecentroid (apply centroid
                           (map geometry/get-vertex
                                (reduce concat edge-points)))
          fcentroid (apply centroid
                           (map geometry/get-vertex
                                face))]
      ($= 4/2 * ecentroid - 2/2 * fcentroid))))

(defn planar-subd-compute-plane-orientation-for-face [struct face]
  (unit-vec
   (apply centroid
          (for [[a b c] (partition 3 (interleave face
                                                 (drop 1 (cycle face))
                                                 (drop 2 (cycle face))))]
            (let [a* (get (:ev->v struct) `(~(edge a b) ~b))
                  c* (get (:ev->v struct) `(~(edge b c) ~b))
                  [A* B C*] (map geometry/get-vertex [a* b c*])]
              (unit-vec ($= (skew-mat ($= C* - B)) <*> ( A* - B))))))))

(defn planar-subd-compute-E [struct face]
  (let [orient (planar-subd-compute-plane-orientation-for-face struct face)
        point  (planar-subd-compute-plane-point-for-face struct face)]
    ($= ( -1 / ( (trans orient) <*> point ) ) * orient)))

(defn planar-subd-compute-E-2 [struct face]
  (let [[a b c d] face
        a- (get (:ev->v struct) `(~(edge a b) ~a))
        a+ (get (:ev->v struct) `(~(edge a b) ~b))
        b- (get (:ev->v struct) `(~(edge b c) ~b))
        b+ (get (:ev->v struct) `(~(edge b c) ~c))
        c- (get (:ev->v struct) `(~(edge c d) ~c))
        c+ (get (:ev->v struct) `(~(edge c d) ~d))
        d- (get (:ev->v struct) `(~(edge d a) ~d))
        d+ (get (:ev->v struct) `(~(edge d a) ~a))
        orient (unit-vec (cross-prod
                          ($= (geometry/get-vertex d-) - (geometry/get-vertex d+))
                          ($= (geometry/get-vertex a+) - (geometry/get-vertex a-))))
        point #_(planar-subd-compute-plane-point-for-face struct face)
        (centroid (map geometry/get-vertex [a- a+ b- b+ c- c+ d- d+]))]
    ($= ( -1 / ( (trans orient) <*> point ) ) * orient)))

(defn planar-subd-face-boundary [struct face]
  (reduce concat
          (for [[a b] (face-to-directed-edges face)]
            (let [v->v  #(get (:ev->v struct) (list (edge a b) %))]
              (list a (v->v a) (v->v b))))))

(defn planar-subd-offset [struct face]
  (let [offset (geometry/get-offset)]
    (geometry/reset-offset)
    (let [[p1 p2 p3 p4 p5 p6
           p7 p8 p9 p10 p11 p12] (doall (map geometry/get-vertex (planar-subd-face-boundary struct face)))
           face-centroid (centroid p2 p3  p5 p6  p8 p9  p11 p12)
           face-normal   (unit-vec (cross-prod ($= p12 - p11) ($= p3 - p2)))]
                                        ;(prn face-centroid face-normal)
      ;(geometry/add-vertex ($= face-centroid - face-normal))
      (geometry/set-offset offset)
      ($= face-centroid - face-normal))))

(defn planar-subd-find-peaks [lst]
  (let [area-per-length (partition 2 (interleave lst (map #(* 10 (/ (quad-area %) (polygon-length %))) lst)))]
    (loop [points   (drop 2 area-per-length)
           previous (first area-per-length)
           point    (second area-per-length)
           acc      '()]
      (let [acc (if (and (< (second previous) (second point))
                         (< (second (first points)) (second point)))
                  (cons (first point) acc)
                  acc)]
        (if (empty? (rest points))
          (do (println (format "Fount %d peaks." (count acc)))
            acc)
          (recur (rest points) point (first points) acc))))))

;; function takes the boundary (vertices)
;; and a list of configurations (vertices).
;; tries to determine the best solution.
(defn select-fixed-point [boundary configs]
  (let [b-length (polygon-length boundary)
        b-area   (quad-area (map #(nth boundary %) '(0 3 6 9)))
        smaller-length  (fn [conf]
                          (< (polygon-length conf) (* 3/4 b-length)))
        smaller-area    (fn [conf]
                          (< (quad-area conf) b-area))]
    ;(prn "Convexity: " (map is-convex? configs))
    ;(prn "Boundary Length: " b-length)
    ;(prn "Boundary Area: " b-area)
    ;(prn "Lengths: " (map polygon-length configs))
    ;(prn "Areas: "(map quad-area configs))
                                        ;(prn "Smaller length: " (map smaller-length configs))

    ;; two cases can occur. (count configs) <= 2 if two distinct fixed
    ;; points are given. or there are more than 2 fixed
    ;; points. E.g. there is a range!
    (println (format ".- FP SELECTOR - Fixed Points: %d -." (count configs)))
    (if (<= (count configs) 2)
      (->> configs
          ;(filter is-convex?)
           (filter smaller-length)
           (filter smaller-area)
           (sort-by quad-area)
           (last))
      (->> configs
           (planar-subd-find-peaks)
           (filter is-convex?)
           (filter smaller-length)
           (filter smaller-area)
           (first)))))

;; comptues a backup point in case of parallelity
;; tries to find the intersection of the line
;; through b and the centroid of (a c) with E.
;; if that fails, checks if all three points are in E;
;; if so it will compute the missing rhombus point.
(defn planar-subd-compute-backup [a b c E]
  (or (line-plane-intersection
       (line-through-points
        b (centroid a c)) E)
      (and (approx= ($= (trans E) <*> a) -1.0 1E-16)
           (approx= ($= (trans E) <*> b) -1.0 1E-16)
           (approx= ($= (trans E) <*> c) -1.0 1E-16)
           ($= a + c - b))))

(defn planar-subd-compute-backup-range-pure [boundary E]
  (let [[a b c] (map #(nth boundary %) [11 0 1])
        offset (planar-subd-compute-backup a b c E)
        length (norm ($= c - a))
        dir    (unit-vec (let [p (plane a b c)]
                                    (if (parrallel? p E)
                                      ($= c - a)
                                      ($= (M (plane a b c)) <*> E))))]
    (for [step (map #(* 1/25 %) (range -50 51))]
      ($= offset + ( step * dir )))))

(defn planar-subd-compute-backup-range [boundary E solve-fn]
  (let [quality (fn [x] (list (is-convex? x) (polygon-length x) (quad-area x)))]
    (for [step (planar-subd-compute-backup-range-pure boundary E)]
      (-> step
          (solve-fn)
          (quality)))))

(defn planar-solution? [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] [q1 q2 q3 q4]]
  (let [ret (every? planar? [[p1 p2 q1 p12]
                             [p2 p3 q2 q1]
                             [p3 p4 p5 q2]
                             [p12 q1 q4 p11]
                             [q1 q2 q3 q4]
                             [q2 p5 p6 q3]
                             [p11 q4 p9 p10]
                             [q4 q3 p8 p9]
                             [p6 p7 p8 q3]])]
    (if-not ret
      (prn (map (juxt planar-dist planar?) [[p1 p2 q1 p12]
                                        [p2 p3 q2 q1]
                                        [p3 p4 p5 q2]
                                        [p12 q1 q4 p11]
                                        [q1 q2 q3 q4]
                                        [q2 p5 p6 q3]
                                        [p11 q4 p9 p10]
                                        [q4 q3 p8 p9]
                                        [p6 p7 p8 q3]])))
    ret))

(defn planar-subd-add-points-and-faces [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] points]
  ;; add new vertices
  (let [[q1 q2 q3 q4] (map #(apply add-vertex %) points)]
    (if (planar-solution? (map get-vertex [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12]) points)
      (->> (for [[a b c d] [[p1 p2 q1 p12]
                            [p2 p3 q2 q1]
                            [p3 p4 p5 q2]
                            [p12 q1 q4 p11]
                            [q1 q2 q3 q4]
                            [q2 p5 p6 q3]
                            [p11 q4 p9 p10]
                            [q4 q3 p8 p9]
                            [q3 p6 p7 p8]]]
             (geometry/add-face a b c d))
           (map face-to-edges)
           (apply concat)
           (set)
           (map (partial apply geometry/add-edge))
           (doall))
      (println "[ E R R O R ] Solution does not generate planar faces!"))))


(defn planar-subd-compute-qs [tp E q1]
  (cond 
   (nil? q1) (prn "Invalid q1")
   (nil? E) (prn "Invalid E")
   (not (and (:a tp) (:b tp) (:c tp) (:d tp) (:f tp) (:g tp) (:h tp) (:i tp))) (prn "Not all lines are available!")
   :else
   (let [q2 (planar-dual/project-onto-plane ($= (M (:c tp)) <*> (M (:b tp)) <*> q1) E)
         q3 (planar-dual/project-onto-plane ($= (M (:i tp)) <*> (M (:f tp)) <*> q2) E)
         q4 (planar-dual/project-onto-plane ($= (M (:g tp)) <*> (M (:h tp)) <*> q3) E)
         q5 (planar-dual/project-onto-plane ($= (M (:a tp)) <*> (M (:d tp)) <*> q4) E)]
     (if-not (approx= 0.0 (norm ($= q5 - q1)) 1E-13)
       (println (format "Diff. of q1 to q5: %2.4e" (norm ($= q5 - q1)))))
     [q2 q3 q4 q5])))


(defn render-geom-highlight-face [boundary face]
  (let [edges (face-to-edges boundary)]
    (doseq [edge edges]
      (remove-edge edge)
      (add-edge (with-meta edge {:color 'red})))
    (with-out-writer (str "/Users/angerman/Dropbox/DA/geom-face-" (join "-" face) ".tikz")
      (println (t/header))
      (doseq [line (render-geometry geometry/*geometry* [5 6 4] [0 0 1] 1.0 13.0 :magnify 15.0 :faces false :edges true :points false)]
        (println line))
      (println (t/footer)))
    (println (format "Wrote geom-face-%s.tikz" (join "-" face)))
    (doseq [edge edges]
      (remove-edge edge)
      (add-edge edge))))

(defn planar-subd-compute-face-points [struct face]
  (let [c (apply centroid (map geometry/get-vertex (planar-subd-face-boundary struct face)))
        dir (planar-subd-compute-plane-orientation-for-face struct face)]
    ;(prn ($= c + (2 * dir)))
    (geometry/set-offset (planar-subd-offset struct face))
    (let [E (planar-subd-compute-E struct face)
          E2 (planar-subd-compute-E-2 struct face)
          boundary (planar-subd-face-boundary struct face)
          vboundary (map geometry/get-vertex boundary)]
      (render-geom-highlight-face boundary face)
      (with-face face
        (println (format "Boundary: %s" (join "-" boundary)))
        #_(prn (apply bind-columns vboundary))
        (prn 'Offset (geometry/get-offset) 'E E)
        (if-let [sol (or (solve vboundary E)
                         (println "---- Second Choice! ----")
                         (solve vboundary E2))]
          (if (planar-solution? vboundary sol)
            (planar-subd-add-points-and-faces boundary sol)
            (println "\n N O T   P L A N A R ! ! ! \n"))
          (println "\n N O   S O L U T I O N ! ! ! \n"))))
    (geometry/reset-offset)
    ))


(defn planar-subd []
  (let [ps (planar-subd-edge-points)
        faces (geometry/get-faces)]
    (def *ps ps)
    (def *faces faces)
    (planar-subd-adjust-vertices ps)
    (let [center (apply centroid (vals (geometry/get-vertices)))]
      (doseq [face faces]
        (geometry/set-offset center)
        (planar-subd-compute-face-points ps face)
        (geometry/reset-offset)
        (apply geometry/remove-face face)
        (doseq [[a b] (face-to-edges face)]
          (geometry/remove-edge a b))))))

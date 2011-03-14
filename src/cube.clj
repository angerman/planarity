(ns cube
  (:use [geometry :only (add-vertex add-face into-faceset-factory parrallel?)])
  (:use [incanter.core :only ($= bind-columns trans)])
  (:require jreality)
  (:require planar-dual)
  (:use [planar-dual :only (M)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)])
  (:use [geometry.shapes.cube])
  (:use [clojure.contrib.duck-streams :only (with-out-writer)]))

(def ifsf (if (bound? #'ifsf)
            ifsf
            (doto (de.jreality.geometry.IndexedFaceSetFactory.)
              (.setGenerateFaceNormals true))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; flags
(def *reset-geom* true)
(def *generate-nodes* true)
(def *use-backup* true)
(def *intersections* false)
(def *solution* 0)


;; legacy
(def *test-e*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Helper Fn

(defn update-geom []
  (into-faceset-factory ifsf))

(defn -main []
  (create-cube)
  (update-geom)
  (jreality/show (jreality/sgc ifsf)))

(defn reset []
  (jreality/dispose)
  (geometry/reset-geom)
  (-main))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Aux Fn


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Planar SubD
;; for each edge, add 1/3, 2/3 points
;; adjust old vertices.

(defn planar-subd-edge-points []
  {:vertices (set (reduce concat (geometry/get-edges)))
   :ev->v (reduce merge
                  (for [[a b] (geometry/get-edges)]
                    (let [A (geometry/get-vertex a)
                          B (geometry/get-vertex b)
                          a* (apply geometry/add-vertex ($= 2/3 * A + 1/3 * B))
                          b* (apply geometry/add-vertex ($= 1/3 * A + 2/3 * B))]
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
        (apply geometry/set-vertex vertex ($= ( 1/2 * vlocation ) + ( 1/2 * ecentroid ))))))) ;; was 1/2, 1/2

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
  (geometry/unit-vec
   (apply centroid
          (for [[a b c] (partition 3 (interleave face
                                                 (drop 1 (cycle face))
                                                 (drop 2 (cycle face))))]
            (let [a* (get (:ev->v struct) `(~(edge a b) ~b))
                  c* (get (:ev->v struct) `(~(edge b c) ~b))
                  [A* B C*] (map geometry/get-vertex [a* b c*])]
              (geometry/unit-vec ($= (geometry/skew-mat ($= C* - B)) <*> ( A* - B))))))))

(defn planar-subd-compute-E [struct face]
  (let [orient (planar-subd-compute-plane-orientation-for-face struct face)
        point  (planar-subd-compute-plane-point-for-face struct face)]
    ($= ( -1 / ( (trans orient) <*> point ) ) * orient)))

(defn planar-subd-face-boundary [struct face]
  (reduce concat
          (for [[a b] (face-to-directed-edges face)]
            (let [v->v  #(get (:ev->v struct) (list (edge a b) %))]
              (list a (v->v a) (v->v b))))))

(defn planar-subd-find-peaks [lst]
  (let [area-per-length (partition 2 (interleave lst (map #(* 10 (/ (quad-area %) (poly-length %))) lst)))]
    (loop [points   (drop 2 area-per-length)
           previous (first area-per-length)
           point    (second area-per-length)
           acc      '()]
      (let [acc (if (and (< (second previous) (second point))
                         (< (second (first points)) (second point)))
                  (cons (first point) acc)
                  acc)]
        (if (empty? (rest points))
          acc
          (recur (rest points) point (first points) acc))))))

;; function takes the boundary (vertices)
;; and a list of configurations (vertices).
;; tries to determine the best solution.
(defn select-fixed-point [boundary configs]
  (let [b-length (poly-length boundary)
        b-area   (quad-area (map #(nth boundary %) '(0 3 6 9)))
        smaller-length  (fn [conf]
                          (< (poly-length conf) (* 3/4 b-length)))
        smaller-area    (fn [conf]
                          (< (quad-area conf) b-area))]
    ;(prn "Convexity: " (map is-convex? configs))
    ;(prn "Boundary Length: " b-length)
    ;(prn "Boundary Area: " b-area)
    ;(prn "Lengths: " (map poly-length configs))
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
  (or (geometry/line-plane-intersection
       (geometry/line-through-points
        b (centroid a c)) E)
      (and (approx= ($= (trans E) <*> a) -1.0 1E-16)
           (approx= ($= (trans E) <*> b) -1.0 1E-16)
           (approx= ($= (trans E) <*> c) -1.0 1E-16)
           ($= a + c - b))))

(defn planar-subd-compute-backup-range-pure [boundary E]
  (let [[a b c] (map #(nth boundary %) [11 0 1])
        offset (planar-subd-compute-backup a b c E)
        length (geometry/norm ($= c - a))
        dir    (geometry/unit-vec (let [p (geometry/plane a b c)]
                                    (if (parrallel? p E)
                                      ($= c - a)
                                      ($= (M (geometry/plane a b c)) <*> E))))]
    (for [step (map #(* 1/25 %) (range -50 51))]
      ($= offset + ( step * dir )))))

(defn planar-subd-compute-backup-range [boundary E solve-fn]
  (let [quality (fn [x] (list (is-convex? x) (poly-length x) (quad-area x)))]
    (for [step (planar-subd-compute-backup-range-pure boundary E)]
      (-> step
          (solve-fn)
          (quality)))))

(def *counter* (atom 0))

(defn planar-subd-add-points-and-faces [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] points]
  ;; add new vertices
  (let [[q1 q2 q3 q4] (map #(apply add-vertex %) points)]
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
         (doall))))

(defn planar-subd-compute-face-points [struct face]
  (let [c (apply centroid (map geometry/get-vertex (planar-subd-face-boundary struct face)))
        dir (planar-subd-compute-plane-orientation-for-face struct face)]
    (geometry/set-offset ($= c - ( 10 * dir )))
    (let [E (planar-subd-compute-E struct face)
          boundary (planar-subd-face-boundary struct face)
          vboundary (map geometry/get-vertex boundary)]
      (if-let [tp (planar-dual/find-first E vboundary)]
        (binding [*test-e* E]
          ;; either choose the backup options range
          ;; if the solution space is 1dimensional
          ;; or choose the fixed points.
          (let [q1-options (or (and *use-backup*
                                    (:equal (meta (:q1 tp)))
                                    (planar-subd-compute-backup-range-pure vboundary E))
                               (filter #($= (trans (:a tp)) <*> %) (:q1 tp)))
                ;; check for `nil` in backups
                ;; backup-options (filter identity [backup-q1])
                quality (fn [v] (poly-length (test-qs tp v)))
                boundary-length (poly-length vboundary)
                configurations (map (partial test-qs tp) q1-options)]
            (if-let [configuration (select-fixed-point vboundary configurations)]
              (planar-subd-add-points-and-faces boundary (cons (last configuration) (butlast configuration)))
              (do ;; else
                (println
                 (apply format "[!! ERROR !!] No configuration found for face %s (EV: %+2.10e; %+2.10e)" (str face) (:evs (meta (:q1 tp)))))
                (if (:equal (meta (:q1 tp)))
                  (with-out-writer (format "/Users/angerman/Dropbox/DA/faces/face-%d.tikz" (swap! *counter* inc))
                    (planar-subd-graph-for-face struct face)))))))
        (prn "ERROR! Could not find first Q1!")))
    (geometry/reset-offset)
    ))





(def *ps)
(def *faces)

(defn planar-subd []
  (let [ps (planar-subd-edge-points)
        faces (geometry/get-faces)]
    (def *ps ps)
    (def *faces faces)
    (planar-subd-adjust-vertices ps)
    (doseq [face faces]
      (planar-subd-compute-face-points ps face)
      (apply geometry/remove-face face)
      (doseq [[a b] (face-to-edges face)]
        (geometry/remove-edge a b)))
    (update-geom)
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Demo Fn
(defn demo-cc []
  (create-cube)
  (update-geom)
  (dotimes [_ 4]
    (sleep)
    (catmull-clark-subd)))

(defn demo-ds []
  (create-cube)
  (update-geom)
  (dotimes [_ 4]
    (sleep)
    (doo-sabin-subd)))

(defn demo-p []
  (create-cube)
  (update-geom)
  (dotimes [_ 2]
    (sleep)
    (planar-subd)
    (sleep)
    (sleep)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; planar stuff.

(def *test-e* '(0 0 -2/3))
(def *test-nodes-ellip*  '((-2 -2 0) (-1 -2 1) (1 -2 1)
                           (2 -2 0) (2 -1 1) (2 1 1)
                           (2 2 0) (1 2 1) (-1 2 1)
                           (-2 2 0) (-2 1 1) (-2 -1 1)))
(def *test-nodes-hyp* '((-2 -2 3) (-1 -2 2) (1 -2 1)
                        (2 -2 0) (2 -1 1) (2 1 2)
                        (2 2 3) (1 2 2) (-1 2 1)
                        (-2 2 0) (-2 1 1) (-2 -1 2)))

(def *test-nodes* *test-nodes-ellip*)

(defn show-nodes []
  (let [nds (map #(apply add-vertex %) *test-nodes*)]
    (doseq [[a b] (face-to-edges nds)]
      (geometry/add-edge a b))
    nds))

(defn show-intersection [[a b c] e]
  (let [h (geometry/plane a b c)
        dir (geometry/unit-vec ($= (geometry/skew-mat h) <*> e))
        pt  (incanter.core/solve (incanter.core/matrix [e h dir]) (incanter.core/matrix [-1 -1 0]))]
    (geometry/add-edge
     (apply add-vertex ($= pt - 2 * dir))
     (apply add-vertex ($= pt + 2 * dir)))))

(defn show-all-intersections [e pts]
  (doseq [d '(2 5 8 11)]
    (show-intersection (take 3 (drop d (cycle pts))) e)))

(defn test-intersections []
  (show-all-intersections *test-e* *test-nodes*))

(defn test-planar []
  (let [E *test-e*
        boundary *test-nodes*]
    (planar-dual/find-first E boundary)))

(defn test-qs [tp q1]
  (let [q2 (planar-dual/project-onto-plane ($= (M (:c tp)) <*> (M (:b tp)) <*> q1) *test-e*)
        q3 (planar-dual/project-onto-plane ($= (M (:i tp)) <*> (M (:f tp)) <*> q2) *test-e*)
        q4 (planar-dual/project-onto-plane ($= (M (:g tp)) <*> (M (:h tp)) <*> q3) *test-e*)
        q5 (planar-dual/project-onto-plane ($= (M (:a tp)) <*> (M (:d tp)) <*> q4) *test-e*)]
    [q2 q3 q4 q5]))

(defn q4->q1 [tp E q4]
  (planar-dual/project-onto-plane ($= (M (:a tp)) <*> (M (:d tp)) <*> q4) E))

(defn q3->q1 [tp E q3]
  (q4->q1 tp E (planar-dual/project-onto-plane ($= (M (:g tp)) <*> (M (:h tp)) <*> q3) E)))

(defn q2->q1 [tp E q2]
  (q3->q1 tp E (planar-dual/project-onto-plane ($= (M (:i tp)) <*> (M (:f tp)) <*> q2) E)))

(declare test-qs)

(defn run-test []
  (if *reset-geom*
    (geometry/reset-geom))
  (if *generate-nodes*
    (do
      (geometry/reset-geom)
      (show-nodes)))
  (let [tp  (test-planar)]
    (if *intersections*
      (test-intersections))
    (let [q1 (if-not (and *use-backup* (:equal (meta (:q1 tp))))
               (nth (filter #($= (trans (:a tp)) <*> %) (:q1 tp)) *solution*)
               ;;use backup point
               (geometry/line-plane-intersection
                (geometry/line-through-points (geometry/get-vertex 0) (centroid (geometry/get-vertex 11)
                                                                                (geometry/get-vertex 1)))
                *test-e*))]
      (if q1
        (let [[q2 q3 q4 q1] (map #(apply add-vertex %) (test-qs tp q1))]
          (geometry/add-edge q1 q2)
          (geometry/add-edge q2 q3)
          (geometry/add-edge q3 q4)
          (geometry/add-edge q4 q1)
          (geometry/add-edge q1 11)
          (geometry/add-edge q1 1)
          (geometry/add-edge q2 2)
          (geometry/add-edge q2 4)
          (geometry/add-edge q3 5)
          (geometry/add-edge q3 7)
          (geometry/add-edge q4 8)
          (geometry/add-edge q4 10))
        (prn "No valid Q1!"))
      (update-geom))))

(defn demo-planar-1 []
  (binding [*use-backup* false *solution* 0] (run-test))
  (sleep)
  (sleep)
  (binding [*use-backup* false *solution* 1] (run-test))
  (sleep)
  (sleep)
  (binding [*use-backup* true] (run-test)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; second planar test
(defn test-planar-2 []
  (geometry/reset-geom)
  (create-cube)
  (let [ps (planar-subd-edge-points)]
    (planar-subd-adjust-vertices ps)
    (def *test-e* (planar-subd-compute-E ps (first (geometry/get-faces))))
    (def *test-nodes* (doall (map geometry/get-vertex (planar-subd-face-boundary ps (first (geometry/get-faces))))))
    (show-nodes)))

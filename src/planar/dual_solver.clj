(ns planar.dual-solver
  (:use [geometry.utils :only (line-plane-intersection
                               line-through-points
                               plane plane-plane-intersection cross-prod
                               impossible? project-onto-plane
                               skew-mat centroid norm unit-vec parrallel? convex?
                               planar-arc
                               quad-area polygon-length avg)])
  (:use [incanter.core :only ($= decomp-eigenvalue abs trans bind-columns to-list)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)])
  (:use [clojure.contrib.seq-utils :only (indexed)])
  (:use [clojure.contrib.duck-streams :only (with-out-writer)])
  (:require tikz)
  (:use [clojure.string :only (join)]))

(def *face*)

(defmacro with-face [face & body]
  `(binding [*face* ~face]
    ~@body))

(defn backup-line [a b c]
  (let [dir ($= c - a)
        off ($= c + a - b)]
    (cross-prod off dir)))

(defn corner-intersection
  "Computes the line of intersection between a corner
   patch (a,b,c) and the hyperplane E.
   If the patch and E are parrallel and non identical,
   returns 'impossible-configuration. Otherwise the line
   through b + (a-b) + (c-b) with the direction of (b-a)"
  [a b c E]
  (let [H (plane a b c)]
    (or (plane-plane-intersection H E)
        (backup-line a b c))))

(defn projection-point
  [a b E]
  (let [off a
        dir ($= b - a)]
    (or (line-plane-intersection [off dir] E)
        dir)))

(defn print-through [args]
  (prn args)
  args)

(defn analyze-eigensystem [es]
  (let [[[val-a vec-a]
         [val-b vec-b]] (->> (interleave (:values es) (trans (:vectors es)))
                             (partition 2)
                             ;; sort by magnitude of
                             ;; eigenvalue.
                             (sort-by (comp abs first))
                             ;; because one eigenvalue is
                             ;; zero. We will drop the
                             ;; smallest!
                             (rest))
         ;; transpose the eigenvectors
         ;; s.t. they are column vectors.
         vec-a (trans vec-a)
         vec-b (trans vec-b)]
    (cond (and (not (approx= val-a 0.0 1E-15))
               (approx= val-a val-b 1E-14))
          (with-meta [vec-a vec-b] {:line true :evs [val-a val-b]})
          (and (not (approx= val-a 0.0 1E-15))
               (not (approx= val-b 0.0 1E-15)))
          (with-meta [vec-a vec-b] {:evs [val-a val-b]})
          (and (not (approx= val-a 0.0 1E-15))
               (approx= val-b 0.0 1E-15))
          (with-meta [vec-a] {:evs [val-a]})
          (and (approx= val-a 0.0 1E-15)
               (not (approx= val-b 0.0 1E-15)))
          (with-meta [vec-b] {:evs [val-b]}))))

(defn pt-mat
  ([lst] (->> lst (map skew-mat) (reduce #($= % <*> %2))))
  ([a & rest] (pt-mat (cons a rest))))

(defn line-solution-offset [a b c E]
  (or (line-plane-intersection
       (line-through-points
        b (centroid a c)) E)
      (and (approx= ($= (trans E) <*> a) -1.0 1E-16)
           (approx= ($= (trans E) <*> b) -1.0 1E-16)
           (approx= ($= (trans E) <*> c) -1.0 1E-16)
           ($= a + c - b))))

(defn line-solution-dir [a b c E]
  (unit-vec (if (parrallel? (plane a b c) E)
              ($= c - a)
              (cross-prod (plane a b c) E))))

(defn line-solution [boundary E]
  (let [[a b c] (map #(nth boundary %) [11 0 1])
        offset (line-solution-offset a b c E)
        length (norm ($= c - a))
        dir    (line-solution-dir a b c E)]
    (for [step (map #(* 1/25 %) (range -50 51))]
      (with-meta 
        ($= offset + (step * dir))
        {:step (* 2.5 step)}))))

(defn solution-quality [solution]
  (* 10 (/ (quad-area solution)
           (polygon-length solution))))

(defn plot-clip []
  (println "\\begin{scope}\\clip (-5,0) rectangle (5,7);"))

(defn plot-end-clip []
  (println "\\end{scope}"))

(defn plot [color data]
  (->> data
       (map (fn [[x y]] (format "(%+2.2f,%+2.3f)" (float x) (float (min y 200)))))
       (join " -- ")
       (format "\\draw[%s] %s;" color)
       (println)))

(defn plot-convex [solutions]
  (print "\\fill[gray!10] ")
  (if (convex? (first solutions))
    (print (format "(%+2.3f,0) " (float (:step (meta (first solutions)))))))
  (doseq [[prev curr] (partition 2 (interleave solutions (rest solutions)))]
    (if (not= (convex? prev)
              (convex? curr))
      (print
       (format
        (if (convex? prev) "rectangle (%+2.3f,7) " "(%+2.3f,0) ")
        (float (avg (map #(:step (meta %)) [prev curr])))))))
  (if (convex? (last solutions))
    (print (format " rectangle (%+2.3f,7)" (float (:step (meta (last solutions)))))))
  (println ";"))

(defn plot-quad-area [solutions]
  (plot "blue" (map (fn [s] [(:step (meta s)) (quad-area s)]) solutions)))

(defn plot-polygon-length [solutions]
  (plot "orange" (map (fn [s] [(:step (meta s)) (polygon-length s)]) solutions)))

(defn plot-quality [solutions]
  (plot "red" (map (fn [s] [(:step (meta s)) (solution-quality s)]) solutions)))

(defn plot-planar-defect [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] solutions]
  (plot "black" (map (fn [s] (let [[q1 q2 q3 q4] s] [(:step (meta s)) (abs (planar-arc p12 q1 q4 p11))])) solutions)))

(declare filter-peaks)

(defn plot-mark-peaks [solutions]
  (->>
   (for [peak (filter-peaks solutions)]
     (format "(%+2.3f,%+2.3f) circle (1pt)"
             (float (:step (meta peak)))
             (float (solution-quality peak))))
   (join " ")
   (format "\\fill[green] %s;")
   (println)))

(defn plot-solutions [boundary b-length b-area solutions]
  (with-out-writer (str "/Users/angerman/Dropbox/DA/face-" (join "-" *face*) ".tikz")
    (println (tikz/header))
    (plot-clip)
    (plot-convex solutions)
    (plot-planar-defect boundary solutions)
    (plot-quad-area solutions)
    (plot-polygon-length solutions)
    (plot-quality solutions)
    (plot-mark-peaks solutions)
    (println (format "\\draw[dashed,orange] (-5,%+2.3f) -- (5,%+2.3f);" b-length b-length))
    (println (format "\\draw[dashed,blue] (-5,%+2.3f) -- (5,%+2.3f);" b-area b-area))
    (plot-end-clip)
    (println (format "\\node at (0,0) [below] {Face: %s};" (join "-" *face*)))
    (println (tikz/footer)))
  (println (format "Wrote face-%s.tikz" (join "-" *face*))))

(defn filter-peaks [range]
;;  (prn "Filteirng Peaks. " (count range))
  (loop [candidates (drop 2 range)
         previous-q (solution-quality (first range))
         candidate  (second range)
         candidate-q nil
         acc        '()]
    (let [current-q (or candidate-q (solution-quality candidate))
          next-q    (solution-quality (first candidates))
          acc       (if (and (< previous-q current-q)
                             (< next-q current-q))
                      (cons candidate acc)
                      acc)]
      (if (empty? (rest candidates))
        acc
        (recur (rest candidates) current-q (first candidates) next-q acc)))))

(defn analyze-planar [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] [q1 q2 q3 q4]]
  (->>
   (map planar-arc  [[p1 p2 q1 p12]
                     [p2 p3 q2 q1]
                     [p3 p4 p5 q2]
                     [p12 q1 q4 p11]
                     [q1 q2 q3 q4]
                     [q2 p5 p6 q3]
                     [p11 q4 p9 p10]
                     [q4 q3 p8 p9]
                     [q3 p6 p7 p8]])
   (map #(format "%+2.2f°" (float %)))
   (join " ")))

(defn planar-defect [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] [q1 q2 q3 q4]]
  (float (planar-arc [p12 q1 q4 p11])))

(defn analyze-planar-only-D [boundary sol]
  (format "%+2.2f°" (planar-defect boundary sol)))

(defn solve
  "Takes the boundary p1,...,p12 and the supporting hyperplane of E,
   and returns q1,q2,q3,q4."
  [boundary E]
;;  (println "Boundary: p1,...,p12")
  ;;  (prn (apply bind-columns E boundary))
  (println "================================================================================")
  (println (format "Solving Face: %s" (join "-" *face*)))
  (let [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] boundary
        a (unit-vec (corner-intersection p12 p1 p2 E))
        b (unit-vec (projection-point    p2 p3 E))
        c (unit-vec (corner-intersection p3 p4 p5 E))
        f (unit-vec (projection-point    p5 p6 E))
        i (unit-vec (corner-intersection p6 p7 p8 E))
        h (unit-vec (projection-point    p8 p9 E))
        g (unit-vec (corner-intersection p9 p10 p11 E))
        d (unit-vec (projection-point    p11 p12 E))]
    (if (not-any? impossible? [a b c d f g h i])
      (let [eigensys (->> [a d g h i f c b]
                          (pt-mat)
                          (decomp-eigenvalue))]
;;        (println "Lines: A, D, G, H, I, F, C, B")
;;        (prn (bind-columns a d g h i f c b))
        (if-let [evs (analyze-eigensystem eigensys)]
          (let [in-E (fn [pt] (project-onto-plane pt E))
                perp-E? (fn [pt] (approx= ($= (trans E) <*> pt) 0.0 1E-16))
                q2 (fn [q1] (in-E ($= (pt-mat c b) <*> q1)))
                q3 (fn [q1] (in-E ($= (pt-mat i f c b) <*> q1)))
                q4 (fn [q1] (in-E ($= (pt-mat g h i f c b) <*> q1)))
                self (fn [q1] ($= (pt-mat a d g h i f c b) <*> q1))
                solution (fn [q1] (with-meta [(in-E q1) (q2 q1) (q3 q1) (q4 q1)]
                                   (meta q1)))
                b-length (polygon-length boundary)
                b-area   (quad-area [p1 p4 p7 p10])
                smaller-length (fn [conf] (< (polygon-length conf)
                                            (* 3/4 b-length)))
                smaller-area   (fn [conf] (< (quad-area conf)
                                            (* 1/1 b-area)))]
            (println (format "Eigenvalues: %s / Diff: %s -- Line: %s"
                             (join "; " (map #(format "%+2.5e" %) (:evs (meta evs))))
                             (if (= 2 (count (:evs (meta evs))))
                               (format "%+2.5e"
                                       (abs (apply - (:evs (meta evs)))))
                               "\\")
                             (str (:line (meta evs)))))
            ;(prn (unit-vec ($= (in-E (second evs)) - (in-E (first evs)))))
            ;(prn (line-solution-dir (last boundary) (first boundary)
                                        ;(secondboundary) E))
            (let [solutions (concat (if (:line (meta evs))
                                      (let [solutions (map solution (line-solution boundary E))]
                                        (plot-solutions boundary b-length b-area solutions)
                                        (filter-peaks solutions)))
                                    (map solution (filter (complement perp-E?) evs)))
                  yes-no-str #(if % "Yes" "No")]
              (println "Sol. | Convex | < Length | < Area |   Length |     Area | Score | Planar Defect | PD ")
              (doseq [[i s] (indexed solutions)]
                (println (format "%4d | %6s | %8s | %6s | %7.2f%% | %7.2f%% | %5.2f | %13s | %s" (inc i) (yes-no-str (convex? s))
                                 (yes-no-str (smaller-length s))
                                 (yes-no-str (smaller-area s))
                                 (float (/ (polygon-length s) b-length 1/100))
                                 (float (/ (quad-area s) b-area 1/100))
                                 (float (solution-quality s))
                                 (analyze-planar-only-D boundary s)
                                 (yes-no-str (complement (approx= (planar-defect boundary s) 0.0 1E-12))))))
              (let [planar-defect-filter (fn [s] (approx= (planar-defect boundary s) 0.0 1E-12))]
                (->> solutions
                     ;(filter convex?)
                     (filter smaller-length)
                     (filter smaller-area)
                     (filter planar-defect-filter)
                     (sort-by quad-area)
                     (last)))))
          ;; else
          (do
            (println "ERROR: Eigensystem posesses no solution.")
            (prn (:values eigensys))))))))

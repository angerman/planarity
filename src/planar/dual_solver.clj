(ns planar.dual-solver
  (:use [geometry.utils :only (line-plane-intersection
                               line-through-points
                               plane plane-plane-intersection cross-prod
                               impossible? project-onto-plane
                               skew-mat centroid norm unit-vec parrallel? convex?
                               planar-arc solution-quality
                               quad-area polygon-length avg)])
  (:use [incanter.core :only ($= decomp-eigenvalue abs trans bind-columns to-list trace sqrt)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)])
  (:use [clojure.contrib.seq-utils :only (indexed)])
  (:use [clojure.contrib.duck-streams :only (with-out-writer)])
  (:require tikz)
  (:use [clojure.string :only (join)]))

(def *face*)
(declare filter-peaks)
(require '[planar.plot :as plot])

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

(defn line-solution [boundary E &{:keys [dir offset] :or {dir nil offset nil}}]
  (let [[a b c] (map #(nth boundary %) [11 0 1])
        offset (or offset (line-solution-offset a b c E))
        length (norm ($= c - a))
        dir    (or dir (line-solution-dir a b c E))]
    (for [step (map #(* 1/25 %) (range -100 101))]
      (with-meta 
        ($= offset + (step * dir))
        {:step (* 1.25 step)}))))


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

(defn is-line? [M E boundary]
  (if-let [off (line-solution-offset (last boundary)
                                     (first boundary)
                                     (second boundary) E)]
    (let [dir (unit-vec (line-solution-dir    (last boundary)
                                              (first boundary)
                                              (second boundary) E))
          x off
          x+ ($= off + dir)
          x- ($= off - dir)]
      (and (parrallel? x ($= M <*> x))
           (parrallel? x+ ($= M <*> x+))
           (parrallel? x- ($= M <*> x-))))
    false))



(def *plot* nil)

(defn solve
  "Takes the boundary p1,...,p12 and the supporting hyperplane of E,
   and returns q1,q2,q3,q4."
  [boundary E]
;;  (println "Boundary: p1,...,p12")
  ;;  (prn (apply bind-columns E boundary))
  (println "================================================================================")
  (println (format "Solving Face: %s" (join "-" *face*)))
  (let [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] boundary
        a (corner-intersection p12 p1 p2 E)
        b (projection-point    p2 p3 E)
        c (corner-intersection p3 p4 p5 E)
        f (projection-point    p5 p6 E)
        i (corner-intersection p6 p7 p8 E)
        h (projection-point    p8 p9 E)
        g (corner-intersection p9 p10 p11 E)
        d (projection-point    p11 p12 E)]
    (if (not-any? impossible? [a b c d f g h i])
      (let [[a b c d f g h i] (map unit-vec [a b c d f g h i])
            M (pt-mat [a d g h i f c b])
            eigensys (decomp-eigenvalue M)]
        ;;        (println "Lines: A, D, G, H, I, F, C, B")
        ;;        (prn (bind-columns a d g h i f c b))
        (let [tr (trace M)
              tr2 (trace ($= M <*> M))]
          ;; det(M-lI) = -l^3 + l^2 tr(M) + l * 1/2 [ tr(M^2) -
          ;; tr^2(M)] + det(M)
          ;; =!= 0
          ;; <=> l(-l^2 + l tr(M) + 1/2 [tr(M^2)-tr^2(M)])
          ;; tr(M)/2 +- sqrt(2tr(M^2) - tr(M)^2)/2
          #_(prn "Traces: " tr ", " tr2)
          (prn ($= tr + (sqrt ($= ( 2 * tr2 ) - ( tr * tr ))))
               ($= tr - (sqrt ($= ( 2 * tr2 ) - ( tr * tr ))))))
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
                                            (* 1/1 b-length))) ;; once
                ;; it was 3/4
                smaller-area   (fn [conf] (< (quad-area conf)
                                            (* 1/1 b-area)))]
            #_(println "Proj. Trans. Matrix")
            #_(prn (format "eigen(t(matrix(nrow=3,c(%s))))"
                         (join ", " (apply concat (to-list (pt-mat a d g h i f c b))))))
            (println "Eigenvalues")
            (prn (:values eigensys))
            (println "Normalized Eigenvectors")
            (prn (apply bind-columns b
                        (for [ev (trans (:vectors eigensys))]
                          (unit-vec (trans ev)))))
            #_(doseq [ev evs]
              (println (format "<ev,e>=%+2.5e" ($= (trans (unit-vec ev)) <*> (unit-vec E))) ))
            (println (format "Eigenvalues: %s / Diff: %s -- Line: %s; %s"
                             (join "; " (map #(format "%+2.5e" %) (:evs (meta evs))))
                             (if (= 2 (count (:evs (meta evs))))
                               (format "%+2.5e"
                                       (abs (apply - (:evs (meta evs)))))
                               "\\")
                             (str (:line (meta evs)))
                             (str (is-line? M E boundary))))
            #_(if (:line (meta evs))
              (do
                (println "Solution is a line. First direction of ev1 -> ev2, second A \\cap E")
                (prn (unit-vec ($= (second evs) - (first evs))))
                (prn (line-solution-dir (last boundary) (first boundary) (second boundary) E))))
            (let [solutions (concat (if (is-line? M E boundary)
                                      (let [solutions (map solution (line-solution boundary E
                                                                                   ;:dir (unit-vec ($= (second evs) - (first evs)))
                                                                                   ;:offset
                                                                                   ;(in-E
                                                                                   ;(first
                                                                                   ;evs))
                                                                                   ))]
                                        (if *plot*
                                          (plot/plot-solutions boundary b-length b-area solutions))
                                        (filter-peaks solutions))
                                      (map solution (filter (complement perp-E?) evs))))
                  yes-no-str #(if % "Yes" "No")
                  planar-defect-filter (fn [s] (approx= (planar-defect boundary s) 0.0 1E-12))]
              (let [sol (->> solutions
                                        ;(filter convex?)
                             (filter smaller-length)
                             (filter smaller-area)
                             (filter planar-defect-filter)
                             (sort-by solution-quality)
                             (last))]
                (println "Sol. | Convex | < Length | < Area |   Length |     Area | Score | Planar Defect | PD ")
                (doseq [[i s] (indexed solutions)]
                  (println (format "%2s%2d | %6s | %8s | %6s | %7.2f%% | %7.2f%% | %5.2f | %13s | %s"
                                   (if (= sol s) "*" " ")
                                   (inc i) (yes-no-str (convex? s))
                                   (yes-no-str (smaller-length s))
                                   (yes-no-str (smaller-area s))
                                   (float (/ (polygon-length s) b-length 1/100))
                                   (float (/ (quad-area s) b-area 1/100))
                                   (float (solution-quality s))
                                   (analyze-planar-only-D boundary s)
                                   (yes-no-str (not (approx= (planar-defect boundary s) 0.0 1E-12))))))
                sol)))
          ;; else
          (do
            (println "ERROR: Eigensystem posesses no solution.")
            (prn (:values eigensys))))))))

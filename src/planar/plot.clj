(ns planar.plot
  (:require tikz)
  (:use [clojure.string :only (join)]
        [clojure.contrib.duck-streams :only (with-out-writer)]
        [geometry.utils :only (avg quad-area polygon-length convex? planar-arc solution-quality)]
        [incanter.core :only (abs)]
        [planar.dual-solver :only (*face* filter-peaks)]))

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

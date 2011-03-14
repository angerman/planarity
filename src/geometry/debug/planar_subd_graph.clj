(ns geometry.debug.planar-subd-graph)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Planar SubD Graphing
(defn planar-subd-graph-header []
  (println "\\begin{tikzpicture}")
  (println "
\\coordinate (p1) at (6,4);
\\coordinate (p2) at (7,4);
\\coordinate (p3) at (8,4);
\\coordinate (p4) at (9,4);
\\coordinate (p5) at (9,5);
\\coordinate (p6) at (9,6);
\\coordinate (p7) at (9,7);
\\coordinate (p8) at (8,7);
\\coordinate (p9) at (7,7);
\\coordinate (p10) at (6,7);
\\coordinate (p11) at (6,6);
\\coordinate (p12) at (6,5);

\\coordinate (q1) at (7,5);
\\coordinate (q2) at (8,5);
\\coordinate (q3) at (8,6);
\\coordinate (q4) at (7,6);

\\foreach \\n/\\o in {A/p1,B/p2,C/p3,D/p12,E/q1,F/q2,G/p11,H/q4,I/q3}
  \\coordinate (\\n) at ($(\\o)+(0.5,0.5)$);

\\draw (p1) -- (p2) -- (p3) -- (p4)
       (p12) -- (q1) -- (q2) -- (p5)
       (p11) -- (q4) -- (q3) -- (p6)
       (p10) -- (p9) -- (p8) -- (p7)
       (p1) -- (p12) -- (p11) -- (p10)
       (p2) -- (q1) -- (q4) -- (p9)
       (p3) -- (q2) -- (q3) -- (p8)
       (p4) -- (p5) -- (p6) -- (p7);
\\node at (p1) [below left]  {$p_1$};
\\node at (p2) [below]       {$p_2$};
\\node at (p3) [below]       {$p_3$};
\\node at (p4) [below right] {$p_4$};
\\node at (p5) [right]       {$p_5$};
\\node at (p6) [right]       {$p_6$};
\\node at (p7) [above right] {$p_7$};
\\node at (p8) [above]       {$p_8$};
\\node at (p9) [above]       {$p_9$};
\\node at (p10)[above left]  {$p_{10}$};
\\node at (p11)[left]        {$p_{11}$};
\\node at (p12)[left]        {$p_{12}$};
\\node at (q1) [below left]  {$q_1$};
\\node at (q2) [below right] {$q_2$};
\\node at (q3) [above right] {$q_3$};
\\node at (q4) [above left]  {$q_4$};

\\foreach \\p in {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,q1,q2,q3,q4}
  \\fill (\\p) circle (2pt);

\\foreach \\n in {A,...,I}
  \\node at (\\n) {$\\n$};

\\node at ($(p7)+(1,0)$) [right] {Schematic depiction of};
\\node at ($(p7)+(1,-0.5)$)   [right] {the $3\\times3$ configuration};
\\node at ($(p7)+(1,-1.5)$) [right] {{\\bf Problem} Find};
\\node at ($(p7)+(1,-2)$) [right] {$q_1,\\ldots,q_4$ such that all};
\\node at ($(p7)+(1,-2.5)$) [right] {faces $A,\\ldots,I$ are planar.};
"))

(defn planar-subd-graph-legend []
  (println "\\fill[orange] (6,3) circle (2pt) node[right=10pt] {Solution Length ($\\lvert\\overline{q_1q_2}\\rvert+\\lvert\\overline{q_2q_3}\\rvert+\\lvert\\overline{q_3q_4}\\rvert+\\lvert\\overline{q_4q_1}\\rvert$)};
\\fill[blue]   (6,2.5) circle (2pt) node[right=10pt] {Solution Area ($\\lvert\\left(q_3-q_1\\right)\\times\\left(q_4-q_2\\right)\\rvert$)};
\\fill[red]    (6,2) circle (2pt) node[right=10pt] {Solution Area/Length (10x)};
\\fill[gray]   (6,1.5) circle (2pt) node[right=10pt] {Convex Solution};
\\fill[green]  (6,1) circle (2pt) node[right=10pt] {Local Maxima};
\\node at (6,0.5) [right] {$x$: direction of $A \\cap E$ with length $\\lvert \\overline{p_{12} p_2} \\rvert$};
\\node at (6,-0.25) [right] {\\small\\emph{dashed lines indicate values for the boundary ($p_1$ --- $p_{12}$)}};
\\node at (5,8.5)  {\\Large Planar Solution Space ($q_1$,$q_2$,$q_3$, and $q_4$)};
\\node at (5,8) {\\small in case of infinite fixed points.};
\\node at (0,-1) {Performance w.r.t. the position of $q_1$};"))

(defn planar-subd-graph-axis [lst]
  (println "\\draw (-5,0) -- (5,0);
\\foreach \\i in {-5,-2.5,0,2.5,5}
 \\draw (\\i,-0.1) -- (\\i,0.1);
\\node at (-5,0) [below] {$-2x$};
\\node at (-2.5,0) [below] {$-1x$};
\\node at (0,0) [below] {$E \\cap \\left( p_1 \\vee \\frac{p_{12}+p_2}{2} \\right)$};
\\node at (2.5,0) [below] {$+1x$};
\\node at (5,0) [below] {$+2x$};")
  lst)

(defn planar-subd-graph-footer []
  (print "\\end{tikzpicture}"))

(defn planar-subd-graph-poly-length [lst]
  (print "\\draw[orange] ")
  (print
   (reduce #(str % " -- " %2 )
           (map (fn [[x y]] (format "(%.3f,%.3f)" (float x) (float (min y 50))))
                (partition 2 (interleave (map #(* 1/10 %) (range -50 51)) (map #(nth % 1) lst))))))
  (print ";\n")
  lst)

(defn planar-subd-graph-area [lst]
  (print "\\draw[blue] ")
  (print
   (reduce #(str % " -- " %2 )
           (map (fn [[x y]] (format "(%.3f,%.3f)" (float x) (float (min y 50))))
                (partition 2 (interleave (map #(* 1/10 %) (range -50 51)) (map #(nth % 2) lst))))))
  (print ";\n")
  lst)

(defn planar-subd-graph-area-per-poly-length [lst]
  (print "\\draw[red] ")
  (print
   (reduce #(str % " -- " %2)
           (map (fn [[x y]] (format "(%.3f,%.3f)" (float x) (float (min y 50))))
                (partition 2 (interleave (map #(* 1/10 %) (range -50 51)) (map #(* 10 (/ (nth % 2) (nth % 1))) lst))))))
  (print ";\n")
  lst)

(defn planar-subd-graph-convex [lst]
  (let [pts (partition 2 (interleave (map #(* 1/10 %) (range -50 51)) (map first lst)))]
    (print "\\fill[gray!10] ")
    (if (second (first pts)) (print (format "(%.3f,0) " (float (first (first pts))))))
    (loop [points   (drop 2 pts)
           point    (second pts)
           previous (first pts)]
      (if (not= (second point)
                (second previous))
        (print
         (if (second previous)
           (format "rectangle (%.3f,7)" (float (avg (first previous) (first point))))
           (format "(%.3f,0)" (float (avg (first previous) (first point)))))))
      (if-not (empty? points)
        (recur (rest points) (first points) point)))
    (if (second (last pts)) (print (format " rectangle (%.3f,7)" (float (first (last pts))))))
    (print ";\n")
    lst))

(defn planar-subd-graph-begin-clip [lst]
  (println "\\begin{scope}")
  (println "\\clip (-5,0) rectangle (5,7);")
  lst)

(defn planar-subd-graph-end-clip [lst]
  (println "\\end{scope}")
  lst)

(defn planar-subd-graph-mark-peaks [lst]
  (print "\\fill[green] ")
  (let [area-per-length (partition 2 (interleave (map #(* 1/10 %) (range -50 51)) (map #(* 10 (/ (nth % 2) (nth % 1))) lst)))]
    (loop [points   (drop 2 area-per-length)
           previous (first area-per-length)
           point    (second area-per-length)]
      (if (and (< (second previous) (second point))
               (< (second (first points)) (second point)))
        (print (apply format "(%.3f,%.3f) circle (2pt) " (map float point))))
      (if-not (empty? (rest points))
        (recur (rest points) point (first points)))))
  (println ";")
  lst)
;; usage (e.g. to render performance for face 10)
;; (def *ps (planar-subd-edge-points))
;; (planar-subd-adjust-vertices *ps)
;; (with-out-writer
;; "/Users/angerman/Dropbox/DA/tikz-render.solution-quality.tikz"
;; (planar-subd-graph-for-face *ps (nth (seq @geometry/*faces*) 10)))

(declare planar-subd-compute-E)
(declare planar-subd-face-boundary)
(declare planar-subd-compute-backup-range)
(declare poly-length)
(declare quad-area)
(declare test-qs)

(defn planar-subd-graph-for-face [ps face]
  (binding [*test-e* (planar-subd-compute-E ps face)]
    (let [boundary (planar-subd-face-boundary ps face)
          vboundary (map geometry/get-vertex boundary)
          tp (planar-dual/find-first *test-e* vboundary)
          b-length (poly-length vboundary)
          b-area   (quad-area (map #(nth vboundary %) '(0 3 6 9)))]
      (planar-subd-graph-header)
      (-> (planar-subd-compute-backup-range vboundary *test-e* (partial test-qs tp))
          (planar-subd-graph-convex)
          (planar-subd-graph-axis)
          (planar-subd-graph-begin-clip)
          (planar-subd-graph-poly-length)
          (planar-subd-graph-area)
          (planar-subd-graph-area-per-poly-length)
          (planar-subd-graph-mark-peaks)
          (planar-subd-graph-end-clip))
      (println (format "\\draw[orange,dashed] (-5,%.3f) -- (5,%.3f);" (float b-length) (float b-length)))
      (println (format "\\draw[blue,dashed] (-5,%.3f) -- (5,%.3f);" (float b-area) (float b-area)))
      (planar-subd-graph-legend)
      (planar-subd-graph-footer))))

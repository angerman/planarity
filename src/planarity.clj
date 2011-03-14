(ns planarity
  (:use [incanter.core :only ($= solve decomp-eigenvalue abs trans matrix bind-rows bind-columns to-list)])
  (:use [geometry.utils :only (plane point norm plane-plane-intersection line-plane-intersection line-through-points skew-mat)])
  (:use [mesh :only (*mesh-ref*)])
  (:use [clojure.contrib.seq-utils :only (indexed)]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The planarity logic. If the tip of the control unit is moved,
;; the tip is interpreted as the dual's hyperplane E equation.
;;   p10 p9  p8  p7
;;    +---+---+---+
;;    | G | H | I |
;;p11 +---+---+---+ p 6
;;    | D | E | F |
;;p12 +---+---+---+ p5
;;    | A | B | C |
;;    +---+---+---+
;;   p1  p2  p3  p4
(def p1  [0 0])
(def p2  [1 0])
(def p3  [2 0])
(def p4  [3 0])
(def p5  [3 1])
(def p6  [3 2])
(def p7  [3 3])
(def p8  [2 3])
(def p9  [1 3])
(def p10 [0 3])
(def p11 [0 2])
(def p12 [0 1])

(defmacro plane-of [mesh a b c]
  `(try
    (apply plane (map #(get-in @~mesh %) [~a ~b ~c]))
    (catch java.lang.IllegalArgumentException e#
      (println (str "Cannot compute dual of " '~a ", " '~b " and " '~c))
      (println (.getMessage e#)))))


;;   p10 p9  p8  p7
;;    +---+---+---+
;;    | G | H | I |
;;p11 +---+---+---+ p 6
;;    | D | E | F |
;;p12 +---+---+---+ p5
;;    | A | B | C |
;;    +---+---+---+
;;   p1  p2  p3  p4

(defn -intersect-hh [[a b c] e]
  "Computes the intersection of a hyperplane through a,b and c and the hyperplane e.
   If both are parrallel no solution is possible and nil will be returned.
   Otherwise a line through b + (a-b) + (c-b) will be returned."
  (let [h (plane a b c)
        l (plane-plane-intersection h e)]
    (cond
     l l
     (and
      (< (abs ($= (trans a) <*> e)) 1e-16)
      (< (abs ($= (trans b) <*> e)) 1e-16)
      (< (abs ($= (trans c) <*> e)) 1e-16)) (plane-plane-intersection
                                  (plane ($= a + c - b)
                                         ($= 2 * c - b)
                                         ($= a + c - b + e))
                                  e)
     :else nil)))

(defn -intersect-lh [[a b] e]
  "Computes the line through a and b intersecting in e.
   If the line is parrallel to e the intersection is equivalent to the pure vector."
  (let [d ($= b - a)]
    (if (< (abs ($= (trans d) <*> e)) 1e-16)
      d
      (line-plane-intersection
       (line-through-points a b)
       e))))

(defn planar-dual [E [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12]]
  (let [a (-intersect-hh [p12 p1 p2] E)
        b (-intersect-lh [p2 p3] E)
        c (-intersect-hh [p3 p4 p5] E)
        f (-intersect-lh [p5 p6] E)
        i (-intersect-hh [p6 p7 p8] E)
        h (-intersect-lh [p8 p9] E)
        g (-intersect-hh [p9 p10 p11] E)
        d (-intersect-lh [p11 p12] E)]
    (prn (bind-columns a b c d f g h i))
    (if (and a b c d f g h i)
      (let [Ma (skew-mat a) Mb (skew-mat b) Mc (skew-mat c) Md (skew-mat d)
            Mf (skew-mat f) Mg (skew-mat g) Mh (skew-mat h) Mi (skew-mat i)
            M (reduce #($= % <*> %2) [Ma Md Mg Mh Mi Mf Mc Mb])
            eig (decomp-eigenvalue M)
            val (map abs (:values eig))
            min-val (apply min val)
            evs (for [[k v] (indexed val) :when (not= v min-val)] k)
            q0  (trans (nth (trans (:vectors eig)) (first evs)))
            q1  ($= Mc <*> Mb <*> q0)
            q2  ($= Mi <*> Mf <*> q1)
            q3  ($= Mg <*> Mh <*> q2)
            [q0 q1 q2 q3] (map #($= % / (-1 * (trans E) <*> %)) [q0 q1 q2 q3])]
        (if (second evs)
          (let [q4  (trans (nth (trans (:vectors eig)) (second evs)))
                q5  ($= Mc <*> Mb <*> q4)
                q6  ($= Mi <*> Mf <*> q5)
                q7  ($= Mg <*> Mh <*> q6)
                [q4 q5 q6 q7] (map #($= % / (-1 * (trans E) <*> %)) [q0 q1 q2 q3])]
            (list [q0 q1 q2 q3]
                  [q4 q5 q6 q7]))
          (list [q0 q1 q2 q3]))))))

(defn make-dual-mesh-planar [mesh E]
  (try
    (let [mesh-at (fn [a b] (get-in @mesh [a b]))
          update-mesh (fn [a b v] (dosync (alter mesh assoc-in [a b] v)))
          a (if-let [A (plane-of mesh p12 p1 p2)]
              (plane-plane-intersection A E))
          b (if-let [B (line-through-points (mesh-at 1 0) (mesh-at 2 0))]
              (line-plane-intersection B E))
          c (if-let [C (plane-of mesh p3 p4 p5)]
              (plane-plane-intersection C E))
          f (if-let [F (line-through-points (mesh-at 3 1) (mesh-at 3 2))]
              (line-plane-intersection F E))
          i (if-let [I (plane-of mesh p6 p7 p8)]
              (plane-plane-intersection I E))
          h (if-let [H (line-through-points (mesh-at 2 3) (mesh-at 1 3))]
              (line-plane-intersection H E))        
          g (if-let [G (plane-of mesh p9 p10 p11)]
              (plane-plane-intersection G E))
          d (if-let [D (line-through-points (mesh-at 0 2) (mesh-at 0 1))]
              (line-plane-intersection D E))]
      ;;    (prn (list a b c f i h g d))
      (if (and a b c d f g h i)
        (let [Ma (skew-mat a)
              Mb (skew-mat b)
              Mc (skew-mat c)
              Md (skew-mat d)
              Mf (skew-mat f)
              Mg (skew-mat g)
              Mh (skew-mat h)
              Mi (skew-mat i)
              M (reduce #($= % <*> %2) [Ma Md Mg Mh Mi Mf Mc Mb])
              eig (decomp-eigenvalue M)
              val (map abs (:values eig))
              min-val (apply min val)
              evs (for [[k v] (indexed val) :when (not= v min-val)] k)
              q0  (trans (nth (trans (:vectors eig)) (first evs)))
              q1  ($= Mc <*> Mb <*> q0)
              q2  ($= Mi <*> Mf <*> q1)
              q3  ($= Mg <*> Mh <*> q2)
              q4  (trans (nth (trans (:vectors eig)) (second evs)))
              q5  ($= Mc <*> Mb <*> q4)
              q6  ($= Mi <*> Mf <*> q5)
              q7  ($= Mg <*> Mh <*> q6)
              [q0 q1 q2 q3
               q4 q5 q6 q7] (map #($= % / (-1 * (trans E) <*> %)) [q0 q1 q2 q3
                                                                   q4 q5 q6 q7])]
          
          
          #_(apply prn (list "A" "C" "E" "G" "I"))
          #_(doseq [point [q0 q1 q2 q3 q4 q5 q6 q7]]
              (apply prn (map #(if (< ($= (trans point) <*> % + 1) 1E-10) "*" " ") [A C E G I]) ))
          (update-mesh 1 1 (to-list q4))
          (update-mesh 2 1 (to-list q5))
          (update-mesh 2 2 (to-list q6))
          (update-mesh 1 2 (to-list q7)))
        (println "Could not compute planer solution."))
      ;;    (update-solution-one q0 q1 q2 q3)
      ;;    (update-solution-two q4 q5 q6 q7)
      ;;    (update-e-mesh b f h d)
      )
    (catch Exception e
      (prn "Failed!"))))

(defn make-primal-mesh-planar [mesh u v]
  (try
    (binding [*mesh-ref* mesh]
      (let [E (val u v)
            rel-val (fn [du dv]
                      (val (+ u du) (+ v dv)))
            rel-set (fn [du dv val]
                      (set (+ u du) (+ v dv) val))]
        (prn (list u v))
        (let [a (plane-plane-intersection (rel-val -1 -1) E)
              b (point (rel-val 0 -2) (rel-val 0 -1) E)
              c (plane-plane-intersection (rel-val  1 -1) E)
              f (point (rel-val 2 0) (rel-val 1 0) E)
              i (plane-plane-intersection (rel-val 1 1) E)
              h (point (rel-val 0 2) (rel-val 0 1) E)
              g (plane-plane-intersection (rel-val -1 1) E)
              d (point (rel-val -2 0) (rel-val -1 0) E)]
          (if (and a b c f i h g d)
            (let [Ma (skew-mat a) Mb (skew-mat b)
                  Mc (skew-mat c) Md (skew-mat d)
                  Mf (skew-mat f) Mg (skew-mat g)
                  Mh (skew-mat h) Mi (skew-mat i)
                  M (reduce #($= % <*> %2) [Ma Md Mg Mh Mi Mf Mc Mb])
                  eig (decomp-eigenvalue M)
                  val (map abs (:values eig))
                  min-val (apply min val)
                  evs (for [[k v] (indexed val) :when (not= v min-val)] k)
                  q0 (trans (nth (trans (:vectors eig)) (first evs)))
                  q1 ($= Mc <*> Mb <*> q0)
                  q2 ($= Mi <*> Mf <*> q1)
                  q3 ($= Mg <*> Mh <*> q2)
                  q4 (trans (nth (trans (:vectors eig)) (second evs)))
                  q5 ($= Mc <*> Mb <*> q4)
                  q6 ($= Mi <*> Mf <*> q5)
                  q7 ($= Mg <*> Mh <*> q6)
                  ;; making sure the points are in E
                  [q0 q1 q2 q3
                   q4 q5 q6 q7] (map #($= % / (-1 * (trans E) <*> %))
                                     [q0 q1 q2 q3 q4 q5 q6 q7])]
              ;; compute the sourrounding planes
              (let [plane-of (fn [l]
                               (apply plane (map #(apply rel-val %) l)))
                    beta    (plane-of [[-1 -1] [-1 -2] [ 0 -2]])
                    gamma   (plane-of [[ 0 -2] [ 1 -2] [ 1 -1]])
                    theta   (plane-of [[ 1 -1] [ 2 -1] [ 2  0]])
                    mu      (plane-of [[ 2  0] [ 2  1] [ 1  1]])
                    omicron (plane-of [[ 1  1] [ 1  2] [ 0  2]])
                    xi      (plane-of [[ 0  2] [-1  2] [-1  1]])
                    iota    (plane-of [[-1  1] [-2  1] [-2  0]])
                    epsilon (plane-of [[-2  0] [-2 -1] [-1 -1]])
                    ;; generic 4 row rhs.
                    rhs     (matrix [-1 -1 -1 -1])
                    slv     (fn [a b c d]
                              (solve (apply bind-rows (map to-list [a b c d])) rhs))]
                ;; now the big solving begins.
                ;; for now we use the solution to the first eigenvalue only.
                (let [B0 (rel-val  0 -1)
                      F0 (rel-val  1  0)
                      H0 (rel-val  0  1)
                      D0 (rel-val -1  0)
                      
                      B (slv q0 q1 beta gamma)
                      F (slv q1 q2 theta mu)
                      H (slv q2 q3 omicron xi)
                      D (slv q3 q0 iota epsilon)

                      B2 (slv q4 q5 beta gamma)
                      F2 (slv q5 q6 theta mu)
                      H2 (slv q6 q7 omicron xi)
                      D2 (slv q7 q4 iota epsilon)

                      diff1 (+ (norm ($= B - B0))
                               (norm ($= F - F0))
                               (norm ($= H - H0))
                               (norm ($= D - D0)))
                      diff2 (+ (norm ($= B2 - B0))
                               (norm ($= F2 - F0))
                               (norm ($= H2 - H0))
                               (norm ($= D2 - D0)))
                      ]
                  (println (format "Solution 1: %f" diff1 ))
                  (println (format "Solution 2: %f" diff2 ))
                  (rel-set  0 -1 (to-list (if (< diff1 diff2) B B2)))
                  (rel-set  1  0 (to-list (if (< diff1 diff2) F F2)))
                  (rel-set  0  1 (to-list (if (< diff1 diff2) H H2)))
                  (rel-set -1  0 (to-list (if (< diff1 diff2) D D2))))))
            (println "Could not compute planar solution")))))
    (catch Exception e
      (do
        (println "Failed!")
        (println (.getMessage e))))))

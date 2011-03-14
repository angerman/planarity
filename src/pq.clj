(ns pq
  (:use [jreality])
  (:use [incanter.core])
  (:use [clojure.contrib.seq-utils :only (indexed)])
  (:import [java.awt Color]))

(def mesh
  (ref [[[-3  3    1] [-1  3  1] [ 1  3  1] [ 3  3  1]]
        [[-3  1    1] [-1  1  1] [ 1  1  1] [ 3  1  1]]
        [[-3 -1    1] [-1 -1  1] [ 1 -1  1] [ 3 -1  1]]
        [[-3 -3    1] [-1 -3  1] [ 1 -3  1] [ 3 -3  1]]]))

(def control-mesh
  [[[4 4 0] [4 5 0]]
   [[5 4 0] [5 5 0]]])

(def e-mesh
  [[[-2  2  0] [ 2  2  0]]
   [[-2 -2  0] [ 2 -2  0]]])

(def solution-one @mesh)

(def solution-two @mesh)

(def line
  (ref [[0 0 0] [1 1 1]]))

(def -mesh (quad-mesh @mesh))

(def -control-mesh (quad-mesh control-mesh))

(def -e-mesh (quad-mesh e-mesh))

(def -solution-one (quad-mesh solution-one))

(def -solution-two (quad-mesh solution-two))

(def -line (line-set @line [0 1]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Mesh update logic
(defn update-mesh [u v p]
  (dosync
   (alter mesh assoc-in [u v] p))
  (doto -mesh
    (.setVertexCoordinates (lst-to-double-array3 @mesh))
    (.update)))

(defn update-e-mesh [a b c d]
  (doto -e-mesh
    (.setVertexCoordinates (lst-to-double-array3 [[a b] [c d]]))
    (.update)))

(defn update-solution-one [a b c d]
  (let [sol @mesh
        sol (assoc-in sol [1 1] a)
        sol (assoc-in sol [2 1] b)
        sol (assoc-in sol [2 2] c)
        sol (assoc-in sol [1 2] d)]
    (doto -solution-one
      (.setVertexCoordinates (lst-to-double-array3 sol))
      (.update))))

(defn update-solution-two [a b c d]
  (let [sol @mesh
        sol (assoc-in sol [1 1] a)
        sol (assoc-in sol [2 1] b)
        sol (assoc-in sol [2 2] c)
        sol (assoc-in sol [1 2] d)]
    (doto -solution-two
      (.setVertexCoordinates (lst-to-double-array3 sol))
      (.update))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Helper function to compute the plane vector (dual)
(defn neg [x] (* -1 x))

(defn skew-mat [[x y z]]
  (matrix [[ 0 (neg z)  y]
           [ z  0 (neg x)]
           [(neg y)  x  0]]))

(defn plane [a b c]
  (solve (bind-rows a b c) (matrix [-1 -1 -1])))

(defn point [p q r] ;; dual version of plane
  (plane p q r))

(defn line-through-points [a b]
  (list a ($= b - a)))

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

(defn mesh-at [u v]
  (get-in @mesh [u v]))

(defn dual-mesh-at [u v]
  (try
    (plane (mesh-at u (inc v))
           (mesh-at u v)
           (mesh-at (inc u) v))
    (catch java.lang.IllegalArgumentException e
      (println (format "Cannot compute dual plane at %d:%d" u v))
      (println (.getMessage e))
      nil)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; A square + a line is used for the orientation controler
(defn update-line [u p]
  (dosync
   (alter line assoc-in [u] p))
  (doto -line
    (.setVertexCoordinates (lst-to-double-array2 @line))
    (.update)))

(defn set-line-to-center []
  (update-line 0 ($= (reduce #($= % + %2)
                             (for [u [0 1]
                                   v [0 1]]
                               (get-in control-mesh [u v])))
                     / 4)))
(defn set-line-tip []
  (update-line 1 ($= (get-in @line [0])
                     + (apply plane (for [u [1 2]
                                          v [2 1]
                                          :when (not= 2 u v)]
                                      (mesh-at u v))))))

(defn update-normal-line []
  (set-line-to-center)
  (set-line-tip))

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

(defmacro plane-of [a b c]
  `(try
    (apply plane (map #(apply mesh-at %) [~a ~b ~c]))
    (catch java.lang.IllegalArgumentException e#
      (println (str "Cannot compute dual of " '~a ", " '~b " and " '~c))
      (println (.getMessage e#)))))

(defn make-mesh-planar [E]
  (let [a (if-let [A (plane-of p12 p1 p2)]
            (plane-plane-intersection A E))
        b (if-let [B (line-through-points (mesh-at 1 0) (mesh-at 2 0))]
            (line-plane-intersection B E))
        c (if-let [C (plane-of p3 p4 p5)]
            (plane-plane-intersection C E))
        f (if-let [F (line-through-points (mesh-at 3 1) (mesh-at 3 2))]
            (line-plane-intersection F E))
        i (if-let [I (plane-of p6 p7 p8)]
            (plane-plane-intersection I E))
        h (if-let [H (line-through-points (mesh-at 2 3) (mesh-at 1 3))]
            (line-plane-intersection H E))        
        g (if-let [G (plane-of p9 p10 p11)]
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
        (update-mesh 1 1 (to-list q0))
        (update-mesh 2 1 (to-list q1))
        (update-mesh 2 2 (to-list q2))
        (update-mesh 1 2 (to-list q3)))
      (println "Could not compute planer solution."))
;;    (update-solution-one q0 q1 q2 q3)
;;    (update-solution-two q4 q5 q6 q7)
;;    (update-e-mesh b f h d)
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The drag tool that controls the top of the orientation line
(defn orient-drag [_ i pos]
  (if (= i 1)
    (do
      (update-line i pos)
      (make-mesh-planar pos))))

(defn mesh-boundary-drag [_ i pos]
  (let [v (mod i 4)
        u (/ (- i v) 4)]
    (if-not (or (= u v 1)
                (= u v 2)
                (= [u v] [1 2])
                (= [u v] [2 1]))
      (update-mesh u v pos))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Debug helper
(defn debug-dual-plane []
  (let [nu (count @mesh)
        nv (count (first @mesh))]
    (doseq [u (range (dec nu))
            v (range (dec nv))]
      (println (format "Dual(%d:%d): %s" u v (str (trans (dual-mesh-at u v))))))))


(defn -main []
  (prn "Here we go!")
  (update-normal-line)
  ;;(debug-dual-plane)
  (make-mesh-planar ($= (nth @line 1) - (nth @line 0)))

  (show (let [cmp (sgc nil)]
          (.addChild cmp (add-drag-fn (sgc -mesh) mesh-boundary-drag))
          ;;(.addChild cmp (sgc -control-mesh))
          ;; (.addChild cmp (sgc -solution-one (red-trans-appearance)))
          ;; (.addChild cmp (sgc -solution-two (green-trans-appearance)))
          ;;(.addChild cmp (sgc -e-mesh))
          ;;(.addChild cmp (add-drag-fn (sgc -line) orient-drag))
          cmp))
  (show (let [cmp (sgc nil)]
          (.addChild cmp (sgc -control-mesh))
          (.addChild cmp (add-drag-fn (sgc -line) orient-drag))
          cmp)))

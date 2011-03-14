(ns pq2
  (:use [jreality])
  (:use [incanter.core])
  (:use [geometry])
  (:use [planarity])
  (:use [clojure.contrib.seq-utils :only (indexed)])
  (:import [java.awt Color]))

(def mesh
  (ref (vec (map vec (partition (count (range -2 3))
                                (for [u (range -2 3)
                                      v (range -2 3)]
                                  [u v 1]))))))


(def e-mesh
  [[[-2  2  0] [ 2  2  0]]
   [[-2 -2  0] [ 2 -2  0]]])

(def solution-one @mesh)

(def solution-two @mesh)

(def -mesh (quad-mesh @mesh))


(def -e-mesh (quad-mesh e-mesh))

(def -solution-one (quad-mesh solution-one))

(def -solution-two (quad-mesh solution-two))

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

(defn make-mesh-planar []
  (let [E (mesh-at 2 2)
        a (plane-plane-intersection (mesh-at 1 1) E)
        b (point (mesh-at 2 0) (mesh-at 2 1) E)
        c (plane-plane-intersection (mesh-at 3 1) E)
        f (point (mesh-at 4 2) (mesh-at 3 2) E)
        i (plane-plane-intersection (mesh-at 3 3) E)
        h (point (mesh-at 2 4) (mesh-at 2 3) E)
        g (plane-plane-intersection (mesh-at 1 3) E)
        d (point (mesh-at 0 2) (mesh-at 1 2) E)]
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
        ;; compute the sourrounding planes
        (let [beta    (plane (mesh-at 1 1) (mesh-at 1 0) (mesh-at 2 0))
              gamma   (plane (mesh-at 2 0) (mesh-at 3 0) (mesh-at 3 1))
              theta   (plane (mesh-at 3 1) (mesh-at 4 1) (mesh-at 4 2))
              mu      (plane (mesh-at 4 2) (mesh-at 4 3) (mesh-at 3 3))
              omicron (plane (mesh-at 3 3) (mesh-at 3 4) (mesh-at 2 4))
              xi      (plane (mesh-at 2 4) (mesh-at 1 4) (mesh-at 1 3))
              iota    (plane (mesh-at 1 3) (mesh-at 0 3) (mesh-at 0 2))
              epsilon (plane (mesh-at 0 2) (mesh-at 0 1) (mesh-at 1 1))
              ;; generic 4 row rhs.
              rhs     (matrix [-1 -1 -1 -1])
              slv     (fn [a b c d] (solve (apply bind-rows (map to-list [a b c d])) rhs))]
          ;; now the big solving begins.
          (let [B (slv q0 q1 beta gamma)
                F (slv q1 q2 theta mu)
                H (slv q2 q3 omicron xi)
                D (slv q3 q0 iota epsilon)]
            (update-mesh 2 1 (to-list B))
            (update-mesh 3 2 (to-list F))
            (update-mesh 2 3 (to-list H))
            (update-mesh 1 2 (to-list D))))
        )
      (println "Could not compute planer solution."))
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The drag tool that controls the top of the orientation line

(defn mesh-boundary-drag [_ i pos]
  (let [nu (count @mesh)
        nv (count (first @mesh))
        v (mod i nv)
        u (/ (- i v) nu)]
    ;; can move any point
    #_(if-not (or (= u v 1)
                (= u v 2)
                (= [u v] [1 2])
                (= [u v] [2 1])))
    (update-mesh u v pos)
    (if (= u v 2)
      (make-primal-mesh-planar mesh u v))))

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
  ;;(debug-dual-plane)
  ;;(make-mesh-planar ($= (nth @line 1) - (nth @line 0)))

  (show (let [cmp (sgc nil)]
          (.addChild cmp (add-drag-fn (sgc -mesh) mesh-boundary-drag))
          ;;(.addChild cmp (sgc -control-mesh))
          ;; (.addChild cmp (sgc -solution-one (red-trans-appearance)))
          ;; (.addChild cmp (sgc -solution-two (green-trans-appearance)))
          ;;(.addChild cmp (sgc -e-mesh))
          ;;(.addChild cmp (add-drag-fn (sgc -line) orient-drag))
          cmp)))

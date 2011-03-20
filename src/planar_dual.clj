(ns planar-dual
  (:use [jreality :only (show add-drag-fn sgc quad-mesh)])
  ;(:use [pq4 :only (init-mesh reset-mesh update-mesh simple-drag)])
  (:use [geometry :only (add-vertex set-vertex add-edge into-edgeset-factory)])
  (:use [geometry.utils :only (plane unit-vec skew-mat)])
  (:use [incanter.core :only (abs $= matrix solve trans decomp-eigenvalue)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)])
  (:require [clojure.string :as str]))


(def mesh)
(declare -mesh)
;(defn update [] (update-mesh -mesh mesh))
;(defn reset [] (reset-mesh mesh))



(def *p1*); []  (get-in @mesh [0 0])
(def *p2*); []  (get-in @mesh [1 0])
(def *p3*); []  (get-in @mesh [2 0])
(def *p4*); []  (get-in @mesh [3 0])
(def *p5*); []  (get-in @mesh [3 1])
(def *p6*); []  (get-in @mesh [3 2])
(def *p7*); []  (get-in @mesh [3 3])
(def *p8*); []  (get-in @mesh [2 3])
(def *p9*); []  (get-in @mesh [1 3])
(def *p10*); [] (get-in @mesh [0 3])
(def *p11*); [] (get-in @mesh [0 2])
(def *p12*); [] (get-in @mesh [0 1])

(defn q1 []  (get-in @mesh [1 1]))
(defn q2 []  (get-in @mesh [2 1]))
(defn q3 []  (get-in @mesh [2 2]))
(defn q4 []  (get-in @mesh [1 2]))

(defn A []   (plane *p12* *p1* *p2*))
(defn C []   (plane *p3* *p4* *p5*))
(defn I []   (plane *p6* *p7* *p8*))
(defn G []   (plane *p9* *p10* *p11*))

;;(defn E []   (plane (q4) (q1) (q2))) ;; This ignores p3!
(def *E*)

(defn plane-plane-intersection [A B]
  (let [a (unit-vec A)
        b (unit-vec B)]
    ;; test if either the vectors are the same or point into the
    ;; oposite direction.
    (if (or (< (reduce + (map #(abs (+ % %2)) a b)) 1e-15)
            (< (reduce + (map #(abs (- % %2)) a b)) 1e-15))
      ;; they are parallel!
      (if (< (reduce + (map #(abs (- % %2)) A B)) 1e-15)
        nil ;; they are identical.
        'impossible-configuration) ;; they are parrallel, but non identical!
      (let [dir ($= (skew-mat a) <*> b)
            off (solve (matrix [A B dir])
                       (matrix [-1 -1 0]))]
        ;; the line is represented by the
        ;; hyperplane through the offset point and
        ;; an arbitrary point =!= offset on the line.
        ($= (skew-mat off) <*> (off + dir))))))

(defn plane-line-intersection [[off dir] P]
  (if-not (< (abs ($= (trans dir) <*> P)) 1e-15) ;; if not parrallel to P
    (let [t ($= ( -1 - (trans off) <*> P) / ( (trans dir) <*> P) )]
      ($= off + t * dir))))

(defn backup-line [a b c]
  (let [dir ($= c - a)
        off ($= c + a - b)]
    ($= (skew-mat off) <*> (off + dir))))

(defmacro gsymb [s]
  (symbol (str "*" s "*")))

(defmacro corner-lines [a b c _ s]
  `(defn ~s []
     (or (plane-plane-intersection (~(symbol (str/capitalize s))) *E*)
         (backup-line (gsymb ~a) (gsymb ~b) (gsymb ~c)))))

(defmacro plane-line-intersections [a b _ s]
  `(defn ~s []
     (or (plane-line-intersection [(gsymb ~a) (map - (gsymb ~b) (gsymb ~a))] *E*)
         (map - (gsymb ~b) (gsymb ~a)))))

(corner-lines p12 p1  p2  -> a)
(corner-lines p3  p4  p5  -> c)
(corner-lines p6  p7  p8  -> i)
(corner-lines p9  p10 p11 -> g)

(plane-line-intersections p2 p3 -> b)
(plane-line-intersections p5 p6 -> f)
(plane-line-intersections p8 p9 -> h)
(plane-line-intersections p11 p12 -> d)


(defn M [v]
  (skew-mat (unit-vec v)))

(defn mmul [& lst]
  (reduce #($= % <*> %2) lst))

(defn project-onto-plane [p h]
  (let [factor ($= (trans h) <*> p)]
    (if (> (abs factor) 1e-17)
      ($= ( -1 / factor ) * p)
      ($= ( -1 / ( (trans h) <*> h) ) * h + p ))))

(defn is-in? [p lst]
  (doseq [l lst]
    (prn ($= (trans p) <*> (trans l))))
  lst)

(defn impossible? [x]
  (= 'impossible-configuration x))

(defn find-q1 []
  (and (empty? (filter impossible? [(a) (b) (c) (d) (f) (g) (h) (i)]))
       (let [es   (decomp-eigenvalue (apply mmul (map #(M (%)) [a d g h i f c b]) ))
             E    *E*
             evs  (take 2 (sort-by abs > (:values es)))
             equal (approx= (first evs) (second evs) 1e-8)]
         
;         (println (str "Eigenvalues: " (str/join ";" (map #(format "%+2.10e" %) evs))))
         (with-meta
           (->> (take 2 (sort-by #(abs (first %)) > 
                                 (partition 2 (interleave (:values es)
                                                          (trans (:vectors es))))))
                (map second)
                ;(is-in? (a))
                (map #(project-onto-plane (trans %) E)))
           {:equal equal
            :evs   evs}))))

(defn proj-trans []
  (let [q2 (apply mmul (map #(M (%)) [c b]))
        q3 (apply mmul (map #(M (%)) [i f c b]))
        q4 (apply mmul (map #(M (%)) [g h i f c b]))
        E  *E*]
    {:q2 (fn [q1] (project-onto-plane ($= q2 <*> q1) E))
     :q3 (fn [q1] (project-onto-plane ($= q3 <*> q1) E))
     :q4 (fn [q1] (project-onto-plane ($= q4 <*> q1) E))}))


(defn find-first [E [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12]]
  (binding [*E* E
            *p1* p1   *p2* p2   *p3* p3   *p4*  p4   *p5*  p5   *p6* p6
            *p7* p7   *p8* p8   *p9* p9   *p10* p10  *p11* p11  *p12* p12]
    (if (empty? (filter impossible? [(a) (b) (c) (d) (f) (g) (h) (i)]))
      {:a  (a) :c (c) :g (g) :i (i)
       :b  (b) :d (d) :f (f) :h (h)
       :q1 (find-q1)
       :pt (proj-trans)})))

(defn find-all []
  (let [q2 (apply mmul (map #(M (%)) [c b]))
        q3 (apply mmul (map #(M (%)) [i f c b]))
        q4 (apply mmul (map #(M (%)) [g h i f c b]))]
    (map (fn [q1] (list q1
                       (project-onto-plane ($= q2 <*> q1) *E*)
                       (project-onto-plane ($= q3 <*> q1) *E*)
                       (project-onto-plane ($= q4 <*> q1) *E*)))
         (find-q1))))


(defn new-e []
  (let [A (unit-vec (A))
        C (unit-vec (C))
        I (unit-vec (I))
        G (unit-vec (G))
        e ($= 1/4 * (A + C + I + G))
        q1 (q1)
        q2 (q2)
        q3 (q3)
        q4 (q4)
        e-off ($= 1/4 * ( q1 + q2 + q3 + q4))
        factor ($= (trans e) <*> e-off)]
    ($= ( -1 / factor ) * e)))

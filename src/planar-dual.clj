(ns planar-dual
  (:use [jreality :only (show add-drag-fn sgc quad-mesh)])
  (:use [pq4 :only (init-mesh reset-mesh update-mesh simple-drag)])
  (:use [geometry :only (plane unit-vec skew-mat add-vertex set-vertex add-edge into-edgeset-factory)])
  (:use [incanter.core :only (abs $= matrix solve trans decomp-eigenvalue)]))


(def mesh)
(declare -mesh)
(defn update [] (update-mesh -mesh mesh))
(defn reset [] (reset-mesh mesh))

(defn -main []
  (def *mesh-size* 4)
  (def mesh (init-mesh *mesh-size*))
  (geometry/reset-geom)
  (def q1a (add-vertex 0 0 0))
  (def q1b (add-vertex 0 0 0))
  (def q2a (add-vertex 0 0 0))
  (def q2b (add-vertex 0 0 0))
  (def q3a (add-vertex 0 0 0))
  (def q3b (add-vertex 0 0 0))
  (def q4a (add-vertex 0 0 0))
  (def q4b (add-vertex 0 0 0))
;;  (add-edge q1a q1b)
;;  (add-edge q2a q2b)
;;  (add-edge q3a q3b)
;;  (add-edge q4a q4b)
  (add-edge q1a q2a)
  (add-edge q2a q3a)
  (add-edge q3a q4a)
  (add-edge q4a q1a)
  
  (add-edge q1b q2b)
  (add-edge q2b q3b)
  (add-edge q3b q4b)
  (add-edge q4b q1b)


  
  (let [-mesh (quad-mesh @mesh)]
    (def iesf (de.jreality.geometry.IndexedLineSetFactory.))
    (into-edgeset-factory iesf)
    (def *qf* -mesh)
    (let [update-fn (partial update-mesh -mesh)
          cmp (sgc nil)]
      (.addChild cmp 
                 (-> (sgc -mesh)
                     (add-drag-fn (partial simple-drag mesh update-fn))))
      (.addChild cmp (sgc iesf))
      (show cmp))))

(defn p1 []  (get-in @mesh [0 0]))
(defn p2 []  (get-in @mesh [1 0]))
(defn p3 []  (get-in @mesh [2 0]))
(defn p4 []  (get-in @mesh [3 0]))
(defn p5 []  (get-in @mesh [3 1]))
(defn p6 []  (get-in @mesh [3 2]))
(defn p7 []  (get-in @mesh [3 3]))
(defn p8 []  (get-in @mesh [2 3]))
(defn p9 []  (get-in @mesh [1 3]))
(defn p10 [] (get-in @mesh [0 3]))
(defn p11 [] (get-in @mesh [0 2]))
(defn p12 [] (get-in @mesh [0 1]))
(defn q1 []  (get-in @mesh [1 1]))
(defn q2 []  (get-in @mesh [2 1]))
(defn q3 []  (get-in @mesh [2 2]))
(defn q4 []  (get-in @mesh [1 2]))

(defn A []   (plane (p12) (p1) (p2)))
(defn C []   (plane (p3) (p4) (p5)))
(defn I []   (plane (p6) (p7) (p8)))
(defn G []   (plane (p9) (p10) (p11)))

(defn E []   (plane (q4) (q1) (q2))) ;; This ignores p3!

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

(defn a []
  (or (plane-plane-intersection (A) (E))
      (backup-line (p12) (p1) (p2))))

(defn c []
  (or (plane-plane-intersection (C) (E))
      (backup-line (p3) (p4) (p5))))

(defn i []
  (or (plane-plane-intersection (I) (E))
      (backup-line (p6) (p7) (p8))))

(defn g []
  (or (plane-plane-intersection (G) (E))
      (backup-line (p9) (p10) (p11))))

(defn b []
  (or (plane-line-intersection [(p2) (map - (p3) (p2))] (E))
      (map - (p3) (p2))))

(defn f []
  (or (plane-line-intersection [(p5) (map - (p6) (p5))] (E))
      (map - (p6) (p5))))

(defn h []
  (or (plane-line-intersection [(p8) (map - (p9) (p8))] (E))
      (map - (p9) (p8))))

(defn d []
  (or (plane-line-intersection [(p11) (map - (p12) (p11))] (E))
      (map - (p12) (p11))))

(defn M [v]
  (skew-mat v))

(defn mmul [& lst]
  (reduce #($= % <*> %2) lst))

(defn project-onto-plane [p h]
  (let [factor ($= (trans h) <*> p)]
    ($= ( -1 / factor ) * p)))

(defn find-q1 []
  (let [es (decomp-eigenvalue (apply mmul (map #(M (%)) [a d g h i f c b])))]
    (prn (:values es))
    (->> (filter #(> (abs (first %)) 1e-15)
                (partition 2 (interleave (:values es)
                                         (:vectors es))))
         (map second)
         (map #(project-onto-plane (trans %) (E))))))

(defn find-first [E [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12]]
  (defn E [] E)
  (defn p1 [] p1)
  (defn p2 [] p2)
  (defn p3 [] p3)
  (defn p4 [] p4)
  (defn p5 [] p5)
  (defn p6 [] p6)
  (defn p7 [] p7)
  (defn p8 [] p8)
  (defn p9 [] p9)
  (defn p10 [] p10)
  (defn p11 [] p11)
  (defn p12 [] p12)
  (find-q1))

(defn find-all []
  (let [q2 (apply mmul (map #(M (%)) [c b]))
        q3 (apply mmul (map #(M (%)) [i f c b]))
        q4 (apply mmul (map #(M (%)) [g h i f c b]))]
    (map (fn [q1] (list q1
                       (project-onto-plane ($= q2 <*> q1) (E))
                       (project-onto-plane ($= q3 <*> q1) (E))
                       (project-onto-plane ($= q4 <*> q1) (E))))
         (find-q1))))

(defn update-vertices []
  (doseq [[n coord] (partition 2 (interleave [q1a q2a q3a q4a] (first (find-all))))]
    (apply set-vertex n coord))
  (doseq [[n coord] (partition 2 (interleave [q1b q2b q3b q4b] (second (find-all))))]
    (apply set-vertex n coord))
  (into-edgeset-factory iesf))

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

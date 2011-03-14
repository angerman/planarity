(ns pq4
  (:use [jreality])
  (:use [geometry.utils :only (point-in-plane plane point unit-vec parallel-plane-through-point line-through-points
                                              line-plane-intersection norm)])
  (:use [planarity :only (make-primal-mesh-planar planar-dual)])
  (:use [incanter.core :only ($= trans matrix sqrt abs bind-columns)]))

(def *mesh-size* 6)

(def mesh)
(def *use-cache* false)
(def hyperplane-cache)

(defn init-mesh [size]
  (ref
   (vec
    (map vec
         (partition size (for [u (range size)
                               v (range size)]
                           [u v 1]))))))
(defn init-cache []
  (def hyperplane-cache
    (ref (vec (map vec
                   (partition (dec *mesh-size*) (for [u (range (dec *mesh-size*))
                                                      v (range (dec *mesh-size*))]
                                                  nil)))))))
(defn reset-cache [cache]
  (doseq [u (range (count @cache))
          v (range (count (first @cache)))]
    (dosync (alter hyperplane-cache assoc-in [u v] nil))))


(declare hyperplane)

(defn get [mesh [u v]]
  (if (and (> u 0) (> v 0) (not= u v 1))
    (try
      (point-in-plane (get-in @mesh [u v])
                      (hyperplane mesh (dec u) (dec v)))
      (catch Exception e
        (println "Failed to compute get...")
        (get-in @mesh [u v])))
    (get-in @mesh [u v])))

(defn set [mesh [u v] val]
  (if-not (nil? val)
    (dosync
     (alter mesh assoc-in [u v] val)
     val) ;; return the input val
    (prn "Cannot set nil value!")))

(defn reset-mesh [mesh]
  (doseq [u (range (count @mesh))
          v (range (count (first @mesh)))]
    (set mesh [u v] [u v 1])))

(defn reset-mesh-nil [mesh]
  (doseq [u (range (count @mesh))
          v (range (count (first @mesh)))
          :when (nil? (get-in @mesh [u v]))]
    (set mesh [u v] [u v 1])))

(defn update-mesh [qf mesh]
  (. qf setVertexCoordinates (lst-to-double-array3 @mesh))
  (. qf update))

;;================================================================================
;; Subdivision

(defn subdivide [mesh]
  (let [nu (count @mesh)
        nv (count (first @mesh))
        U  (- (* 3 nu) 2)
        V  (- (* 3 nv) 2)
        new-mesh (ref (vec (map vec (partition V (for [u (range U)
                                                       v (range V)]
                                                   [(/ u 3) (/ v 3) 1])))))
        get (fn [mesh [u v]] (get-in @mesh [u v]))]
    (prn "Copying pts")
    ;; set all older points.
    (doseq [u (range nu)
            v (range nv)]
      (set new-mesh [(* 3 u) (* 3 v)] (get mesh [u v])))
    (prn "Computing intermediate pts")
    ;; compute intermediate points
    (doseq [u (range (- nu 1))
            v (range (- nv 1))]
      (set new-mesh [(+ (* 3 u) 1) (* 3 v)] ($= 2/3 * (get mesh [u v]) + 1/3 * (get mesh [(+ u 1) v])))
      (set new-mesh [(+ (* 3 u) 2) (* 3 v)] ($= 1/3 * (get mesh [u v]) + 2/3 * (get mesh [(+ u 1) v])))
      (set new-mesh [(* 3 u) (+ (* 3 v) 1)] ($= 2/3 * (get mesh [u v]) + 1/3 * (get mesh [u (+ v 1)])))
      (set new-mesh [(* 3 u) (+ (* 3 v) 2)] ($= 1/3 * (get mesh [u v]) + 2/3 * (get mesh [u (+ v 1)])))
      (set new-mesh [(+ (* 3 u) 1) (+ (* 3 v) 1)] ($= 1/3 * (get mesh [u v])
                                                      + 1/3 * (get mesh [(+ u 1) v])
                                                      + 1/3 * (get mesh [u (+ v 1)])))
      (set new-mesh [(+ (* 3 u) 2) (+ (* 3 v) 1)] ($= 1/3 * (get mesh [(+ u 1) v])
                                                      + 1/3 * (get mesh [u v])
                                                      + 1/3 * (get mesh [(+ u 1) (+ v 1)])))
      (set new-mesh [(+ (* 3 u) 2) (+ (* 3 v) 2)] ($= 1/3 * (get mesh [(+ u 1) (+ v 1)])
                                                      + 1/3 * (get mesh [(+ u 1) v])
                                                      + 1/3 * (get mesh [u (+ v 1)])))
      (set new-mesh [(+ (* 3 u) 1) (+ (* 3 v) 2)] ($= 1/3 * (get mesh [u (+ v 1)])
                                                      + 1/3 * (get mesh [(+ u 1) (+ v 1)])
                                                      + 1/3 * (get mesh [u v]))))
    (prn "Computing boundary")
    ;; compute the last border
    (doseq [u (range (- nu 1))]
      (set new-mesh [(+ (* 3 u) 1) (* 3 (- nv 1))] ($= 2/3 * (get mesh [u (- nv 1)]) + 1/3 * (get mesh [(+ u 1) (- nv 1)])))
      (set new-mesh [(+ (* 3 u) 2) (* 3 (- nv 1))] ($= 1/3 * (get mesh [u (- nv 1)]) + 2/3 * (get mesh [(+ u 1) (- nv 1)]))))
    (doseq [v (range (- nv 1))]
      (set new-mesh [(* 3 (- nu 1)) (+ (* 3 v) 1)] ($= 2/3 * (get mesh [(- nu 1) v]) + 1/3 * (get mesh [(- nu 1) (+ v 1)])))
      (set new-mesh [(* 3 (- nu 1)) (+ (* 3 v) 2)] ($= 1/3 * (get mesh [(- nu 1) v]) + 2/3 * (get mesh [(- nu 1) (+ v 1)]))))

    (prn "Updaing pts")
    ;; smoothing
    (doseq [u (range 1 (- nu 1))
            v (range 1 (- nv 1))]
      (set new-mesh [(* 3 u) (* 3 v)] ($= 2/6 * (get mesh [u v])
                                          + 1/6 * (get new-mesh [(- (* 3 u) 1) (* 3 v)])
                                          + 1/6 * (get new-mesh [(* 3 u) (- (* 3 v) 1)])
                                          + 1/6 * (get new-mesh [(+ (* 3 u) 1) (* 3 v)])
                                          + 1/6 * (get new-mesh [(* 3 u) (+ (* 3 v) 1)]))))
    ;; boundary smoothing.
    (doseq [u (range 1 (- nu 1))]
      (let [v (- nv 1)]
        (set new-mesh [(* 3 u) (* 3 v)] ($= 2/4 * (get mesh [u v])
                                            + 1/4 * (get new-mesh [(- (* 3 u) 1) (* 3 v)])
                                            + 1/4 * (get new-mesh [(+ (* 3 u) 1) (* 3 v)])))))
    (doseq [v (range 1 (- nv 1))]
      (let [u (- nu 1)]
        (set new-mesh [(* 3 u) (* 3 v)] ($= 2/4 * (get mesh [u v])
                                            + 1/4 * (get new-mesh [(* 3 u) (- (* 3 v) 1)])
                                            + 1/4 * (get new-mesh [(* 3 u) (+ (* 3 v) 1)])))))
    new-mesh))

;;================================================================================
;; copy-patch copies a 3x3 patch from an existing mesh
(defn copy-path [mesh u v]
  (vec (map vec
            (let [su (* 3 u)
                  sv (* 3 u)]
              (partition 4 (for [u  (range 4)
                                 v  (range 4)]
                             (get-in @mesh [(+ su u) (+ sv v)])))))))

(defn is-planar? [mesh u v]
  (let [h (apply plane (map #(get-in @mesh %) [[u v] [(inc u) v] [u (inc v)]]))]
    (< (abs ($= 1 + (trans h) <*> (get-in @mesh [(inc u) (inc v)]))) 1E-12)))

(defn planarity-info [mesh]
  (reduce #(str % "\n" %2)
          (map #(apply str %)
               (partition 3
                          (for [u (range 3)
                                v (reverse (range 3))]
                            (if (is-planar? mesh u v) "[*]" "[ ]"))))))

;; expects a 3x3 mesh and computes whether or not the algorithm can work.
(defn is-applicable? [mesh]
  (let [E (apply plane (map #(get-in @mesh %) [1 2] [1 1] [2 1]))])
  )
;;================================================================================

(defn plane-of [mesh a b c]
  (apply plane (map #(get mesh %) [a b c])))

(defn hyperplane [mesh u v]
  (if (or (not *use-cache*) (nil? (get-in @hyperplane-cache [u v])))
    (let [h
          (if (or (= u 0) (= v 0))
            (plane-of mesh [u (+ v 1)] [u v] [(+ u 1) v])
            (let [I   (unit-vec (hyperplane mesh (- u 1) v))
                  II  (unit-vec (hyperplane mesh (- u 1) (- v 1)))
                  III (unit-vec (hyperplane mesh u (- v 1)))
                  n   (unit-vec ($= I + III))
                  IV  (unit-vec ($= ( 2 * (trans n) <*> II) * n - II ))]
              (parallel-plane-through-point IV (get mesh [u v]))))]
      (dosync (alter hyperplane-cache assoc-in [u v] h))
      h)
    (get-in @hyperplane-cache [u v])))

(def *qf*)

(defn compute-u [u v]
  (if (< u (- (count @mesh) 1))
    (try
      (point
       (hyperplane mesh (- u 1) v)
       (hyperplane mesh (- u 1) (- v 1))
       (hyperplane mesh u (- v 1)))
      (catch IllegalArgumentException e ;; assuming Matrix is
        ;; singular. -- all hyperplanes are the same
        ($= (get mesh [(- u 1) v])
            + (get mesh [u (- v 1)])
            - (get mesh [(- u 1) (- v 1)]))))
    (let [dir (line-through-points (get mesh [u (- v 1)])
                                   (get mesh [u v]))]
      (line-plane-intersection dir (hyperplane mesh (dec u) v)))))

(defn compute-v [u v]
  (if (< v (- (count (first @mesh)) 1))
    (try
      (point
       (hyperplane mesh (- u 1) v)
       (hyperplane mesh (- u 1) (- v 1))
       (hyperplane mesh u (- v 1)))
      (catch IllegalArgumentException e ;; assuming Matrix is
        ;; singular. -- all hyperplanes are the same
        ($= (get mesh [(- u 1) v])
            + (get mesh [u (- v 1)])
            - (get mesh [(- u 1) (- v 1)])))))
    (let [dir (line-through-points (get mesh [(- u 1) v])
                                   (get mesh [u v]))]
      (line-plane-intersection dir (hyperplane mesh u (dec v)))))

(defn compute-all-u [v]
  (loop [u (inc v)]
    (set mesh [u v] (compute-u u v))
    (if (< u (dec (count @mesh)))
      (recur (inc u)))))

(defn compute-all-v [u]
  (loop [v (inc u)]
    (set mesh [u v] (compute-v u v))
    (if (< v (dec (count (first @mesh))))
      (recur (inc v)))))

(defn blah []
  (let [IV (let [I   (unit-vec (plane-of mesh [0 2] [0 1] [1 1]))
                 II  (unit-vec (plane-of mesh [0 1] [0 0] [1 0]))
                 III (unit-vec (plane-of mesh [1 1] [1 0] [2 0]))
                 n   (unit-vec ($= I + III))
                 IV  (unit-vec ($= ( 2 * (trans n) <*> II) * n - II))
                 scale ($= -1 * (trans IV) <*> (get-in @mesh [1 1]))]
             ($= IV / scale ))]
    (prn IV)
    (dosync
     (alter mesh assoc-in [2 1] (line-plane-intersection
                                 (line-through-points (get-in @mesh [2 0]) (get-in @mesh [2 1]))
                                 IV))
     (alter mesh assoc-in [1 2] (line-plane-intersection
                                 (line-through-points (get-in @mesh [0 2]) (get-in @mesh [1 2]))
                                 IV)))
    nil))


(defn compute-all []
  (binding [*use-cache* true]
    (reset-cache hyperplane-cache)
    ;; first compute first row of u and v
    (compute-all-u 1)
    (compute-all-v 1)
    ;; compute the point one further in (2 2)
    (set mesh [2 2] (point (hyperplane mesh 1 2)
                           (hyperplane mesh 1 1)
                           (hyperplane mesh 2 1))))
  (update-mesh *qf* mesh))

(defn make-voss []
  (try
    (binding [*use-cache* true]
      (reset-cache hyperplane-cache)
      ;; we assume [1 1] is set (e.g. moved by hand)
      (doseq [n (range 2 (dec (count @mesh)))]
        ;; (println (format "Making %d..." n))
        (compute-all-u (dec n))
        (compute-all-v (dec n))
        (set mesh [n n] (point
                         (hyperplane mesh (dec n) n)
                         (hyperplane mesh (dec n) (dec n))
                         (hyperplane mesh n (dec n)))))
      (let [n (dec (count @mesh))]
        (doseq [m (range 2 n)]
          (set mesh [n m] (get mesh [n m]))
          (set mesh [m n] (get mesh [m n])))
        (set mesh [n n] (get mesh [n n]))))
    (catch Exception e
      (do
        (println "Voss Failed!")
        (println (.getMessage e))
        (.printStackTrace e))))
  (update-mesh *qf* mesh))


(defn voss-boundary-move [mesh updatefn _ i pos]
  (let [nu (count @mesh)
        nv (count (first @mesh))
        v (mod i nv)
        u (/ (- i v) nu)]
    (if (or (= u 0) (= v 0)
            (= u (dec nu)) (= v (dec nv)))
      (do
        (set mesh [u v] pos)
        (make-voss)
        ;;(updatefn mesh)
        ))))

(defn voss-angle-move [mesh updatefn _ i pos]
  (let [nu (count @mesh)
        nv (count (first @mesh))
        v (mod i nv)
        u (/ (- i v) nu)
        rel-get (fn [du dv] (get mesh [(+ u du) (+ v dv)]))]
    (if (and
         (or (= u 1) (= v 1))
         (not= v 0) (not= u 0))
      (let [a (rel-get -1 0)
            b (rel-get -1 -1)
            c (rel-get 0 -1)]
        (if (and a b c)
          (let [new-pos (point-in-plane pos (hyperplane mesh (dec u) (dec v)))]
            (do
              (set mesh [u v] new-pos)
              (make-voss)
              #_(if (= u v 1)
                (make-voss)
                (updatefn mesh))
              )))))))

(defn simple-drag [mesh updatefn _ i pos]
  (let [nu (count @mesh)
        nv (count (first @mesh))
        v  (mod i nv)
        u  (/ (- i v) nu)]
    (do
      (set mesh [u v] pos)
      (updatefn mesh))))

(defn knod-drag [mesh updatefn _ i pos]
  (let [nu (count @mesh)
        nv (count (first @mesh))
        v  (mod i nv)
        u  (/ (- i v) nu)]
    (if (and (> u 1) (> v 1) (< u (- nu 2)) (< v (- nv 2)))
      (do
        (set mesh [u v] pos)
        (make-primal-mesh-planar mesh u v)
        (updatefn mesh)))))

(defn planar-fn [mesh updatefn _ i pos]
  (let [p1 (get-in @mesh [0 0])
        p2 (get-in @mesh [1 0])
        p3 (get-in @mesh [2 0])
        p4 (get-in @mesh [3 0])
        p5 (get-in @mesh [3 1])
        p6 (get-in @mesh [3 2])
        p7 (get-in @mesh [3 3])
        p8 (get-in @mesh [2 3])
        p9 (get-in @mesh [1 3])
        p10 (get-in @mesh [0 3])
        p11 (get-in @mesh [0 2])
        p12 (get-in @mesh [0 1])
        o0 (get-in @mesh [1 1])
        o1 (get-in @mesh [2 1])
        o2 (get-in @mesh [2 2])
        o3 (get-in @mesh [1 2])]
    (let [[[q0 q1 q2 q3] & rest] (planar-dual (plane (get-in @mesh [1 2])
                                                            (get-in @mesh [1 1])
                                                            (get-in @mesh [2 1]))
                                              [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12])]
      (prn (bind-columns q0 q1 q2 q3))
      (if rest
        (let [[q4 q5 q6 q7] (first rest)
              _ (prn rest)
              diff1 (reduce + (map norm [($= o0 - q0) ($= o1 - q1) ($= o2 - q2) ($= o3 - q3)]))
              _ (prn diff1)
              diff2 (reduce + (map norm [($= o0 - q4) ($= o1 - q5) ($= o2 - q6) ($= o3 - q7)]))
              _ (prn diff2)]
          (do
            (prn (list diff1 diff2))
            (set mesh [1 1] (if (< diff1 diff2) q0 q4))
            (set mesh [2 1] (if (< diff1 diff2) q1 q5))
            (set mesh [2 2] (if (< diff1 diff2) q2 q6))
            (set mesh [1 2] (if (< diff1 diff2) q3 q7))
            ;;        (updatefn mesh)
            ))
        (do
          (prn "Only one!")
          (set mesh [1 1] q0)
          (set mesh [2 1] q1)
          (set mesh [2 2] q2)
          (set mesh [1 2] q3))))))

(defn reveal-origin []
  (set mesh [0 0] [0 0 2])
  (update-mesh *qf* mesh))


(defn go []
  (let [-mesh (quad-mesh @mesh)]
    (def *qf* -mesh)
    (let [update-fn (partial update-mesh -mesh)]
      (show (-> (sgc -mesh)
                (add-drag-fn (partial voss-angle-move mesh update-fn))
                (add-drag-fn (partial voss-boundary-move mesh update-fn))
                (add-drag-fn (partial knod-drag mesh update-fn)))))))

(defn -main []
  (init-mesh)
  (init-cache)
  (go)
  (add-key-event (last-viewer) :S reveal-origin)
  (add-key-event (last-viewer) :R (partial reset-mesh mesh))
  (add-key-event (last-viewer) :V make-voss))

(defn -test []
  (let [-mesh (quad-mesh @mesh)]
    (def *qf* -mesh)
    (let [update-fn (partial update-mesh -mesh)]
      (show (-> (sgc -mesh)
                (add-drag-fn (partial simple-drag mesh update-fn)))))))

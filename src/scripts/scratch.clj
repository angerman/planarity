;; extracting one of the edge patches from the
;; planar cube subdivision

(require '[geometry :as g])
(require '[jreality :as jr])
(require '[geometry.shapes.cube :as cube])
(require '[geometry.subdivision.planar :as planar])
(def *isfs* (jr/indexed-face-set-factory))
(def *cube* (cube/create))
(jr/show (jr/sgc *isfs*))

;; helper fn
(defn update-geom [] (g/with-geometry *cube* (g/into-faceset-factory *isfs*)))
;; show cube
(update-geom)

(g/with-geometry *cube* (planar/planar-subd))
;; show cube lvl 1
(update-geom)

(g/with-geometry *cube* (planar/planar-subd))
;; show cube lvl 2
(update-geom)

;; edge boundary and corresponding vertices
(def *b* (map #(Integer/parseInt %) (clojure.string/split "32-256-255-5-161-162-13-239-240-34-128-127" #_"32-272-271-5-167-168-13-253-254-35-152-151" #"-")))
(def *vb* (g/with-geometry *cube* (doall (map g/get-vertex *b*))))

;; new geometry, add vertices and edges
(defn reset-boundary []
  (def *boundary* (g/create-geometry))
  (g/with-geometry *boundary* (doseq [v *vb*] (g/add-vertex v)))
  (g/with-geometry *boundary* (let [nds (g/get-nodes)] (doseq [e (partition 2 1 [ (first nds) ] nds)] (g/add-edge e)))))
(reset-boundary)
;; show
(defn update-boundary []
  (g/with-geometry *boundary* (g/into-faceset-factory *isfs*)))

(update-boundary)

(defn shift [n coll]
  (let [l (count coll)]
    (take l (drop (mod n l) (cycle coll)))))

(use '[incanter.core :only ($= trans matrix)])
(use '[geometry.utils :only (centroid unit-vec cross-prod)])

(defn e-dir-one []
  (centroid (g/with-geometry *boundary*
              (for [[a b c] (partition 3 (doall
                                          (map g/get-vertex (shift -1 (g/get-nodes)))))]
                (unit-vec (cross-prod ($= a - b) ($= c - b)))))))

(defn e-dir-two []
  (g/with-geometry *boundary*
    (let [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] (map g/get-vertex (g/get-nodes))]
      (unit-vec (cross-prod ($= p11 - p12) ($= p3 - p2))))))

(def lambda 0.5)
(defn e-dir []
  (map +
       (map #(* lambda %) (e-dir-one))
       (map #(* (- 1 lambda) %) (e-dir-two))))

(defn e-pt []
  (g/with-geometry *boundary*
    (let [V (partition 3 (doall (map g/get-vertex (g/get-nodes))))
          A (centroid (map first V))
          B (centroid (apply concat (map rest V)))]
      [A B])))

(defn e-offset [x]
  (let [[A B] (e-pt)]
    ($= ( x * B ) + ( ( 1 - x ) * A))))

(defn e [x]
  (let [orient (unit-vec (e-dir))
        point  (e-offset x)]
    ($= ( -1 / ( (trans orient) <*> point ) ) * orient)))

(use '[planar.dual-solver :only (solve with-face) ])


;; utils to show the intersections
(defn show-intersection [[a b c] e]
  (let [h (geometry.utils/plane a b c)
        dir (unit-vec (cross-prod h e))
        pt  (geometry.utils/line-plane-intersection [b ($= (centroid a c) - b)] e)]
    (g/add-vertex pt)
    (g/add-edge
     (apply g/add-vertex ($= pt - 1 * dir))
     (apply g/add-vertex ($= pt + 1 * dir)))))

(defn show-all-intersections [e pts]
  (doseq [d '(2 5 8 11)]
    (show-intersection (take 3 (drop d (cycle pts))) e)))





;; original data
;; Boundary: 32-256-255-5-161-162-13-239-240-34-128-127
;; [-1.4568 -1.4259 -1.3519 -1.2531 -1.0185 -0.7593 -0.5123 -0.5741 -0.6481 -0.6728 -0.9815 -1.2407
;; -1.4568 -1.4259 -1.3519 -1.2531 -1.3519 -1.4259 -1.4568 -1.5741 -1.6481 -1.6975 -1.6481 -1.5741
;;  0.5123  0.7593  1.0185  1.2531  1.3519  1.4259  1.4568  1.2407  0.9815  0.6728  0.6481  0.5741]

;; (1.4999999999999991 1.5000000000000002 1.4999999999999998)
;; (0.17056226070096847 0.4236943594242576 -0.17056226070096833)
;; ================================================================================
;; Solving Face: 32-5-13-34
;; Eigenvalues: +1.82880e-01; +1.82880e-01 / Diff: +0.00000e+00 -- Line: true
;; Wrote face-32-5-13-34.tikz
;; Sol. | Convex | < Length | < Area |   Length |     Area | Score | Planar Defect | PD 
;;    1 |    Yes |      Yes |    Yes |   37.80% |   11.90% |  0.59 |        +3.30° | Yes
;;    2 |     No |       No |     No | 23788.32% | 156715.41% | 12.44 |        +0.00° | Yes
;;    3 |     No |       No |     No | 3391427444736.00% | 153529.83% |  0.00 |        -0.00° | No
;; Wrote geom-face-21-56-13-5.tikz

;; set offset and get local boundary. // WARING: DIFFERENT OFFSETS ->
;; Different EIGENVALUES+EIGENVECTORS!! Possibly Complex!
(g/with-geometry *boundary* (g/set-offset 1.5 1.5 1.5))
(defn lvb [] (g/with-geometry *boundary* (doall (for [i (range 1 13)] (g/get-vertex i)))))

;; compute solution
(dotimes [i 21]
  (let [dist (+ 0 (/ (- i 10) 5))]
    (with-face `(0 0 0 ~(float dist)) (solve (lvb) (e dist)))))

;; TODO: What happens if we "reverse" the boundary?
;; Or shift it by 3,6,9?





;;================================================================================
;;================================================================================
;; shape for descriptive geometry
(defn dg-shape []
  (let [geom (g/create-geometry)]
    (g/with-geometry geom
      (let [;; first layer (x = 0)
            a (g/add-vertex 0 1 1)
            b (g/add-vertex 0 0 1)
            c (g/add-vertex 0 0 3)
            d (g/add-vertex 0 2 3)
            e (g/add-vertex 0 2 1)
            f (g/add-vertex 0 4 1)
            g (g/add-vertex 0 4 0)
            h (g/add-vertex 0 1 0)
            ;; second layer (x = 2)
            i (g/add-vertex 2 2 3)
            j (g/add-vertex 2 2 1)
            k (g/add-vertex 2 4 1)
            l (g/add-vertex 2 4 2)
            m (g/add-vertex 2 3 3)
            ;; third layer (x = 3)
            o (g/add-vertex 3 1 1)
            p (g/add-vertex 3 0 1)
            q (g/add-vertex 3 0 3)
            r (g/add-vertex 3 3 3)
            s (g/add-vertex 3 4 2)
            t (g/add-vertex 3 4 0)
            u (g/add-vertex 3 1 0)]
        (doseq [e [[a b] [b c] [c d] [d e] [e f] [f g] [g h] [h a] ;; first layer boundary
                   [i j] [j k] [k l] [l m] [m i]                   ;; second layer boundary
                   [o p] [p q] [q r] [r s] [s t] [t u] [u o]       ;; third layer boundary
                   [a o] [b p] [c q] [g t] [h u]                   ;; outer hull (first <-> third layer)
                   [d i] [e j] [f k]                               ;; inner edges: first <-> second layer
                   [l s] [m r]]]                                   ;; inner edges: second <-> third layer
          (g/add-edge e))))
    geom))

(doseq [l (concat
           (r/render-geometry *dg* [4   2 1.5] [ 0 0 1] 0.5 5 :points false :faces false :parallel true :target [1.5 2 1.5])
           (r/render-geometry *dg* [1.5 2 4  ] [-1 0 0] 0.5 5 :points false :faces false :parallel true :target [1.5 2 1.5] :offset [0 4])
           (r/render-geometry *dg* [1.5 6 1.5] [-1 0 0] 0.5 5 :points false :faces false :parallel true :target [1.5 2 1.5] :offset [4 4]))]
  (println l))

(doseq [l (r/render-geometry *dg* [-2 5 5] [0 0 1] 0.5 10 :points false :target [1.5 2 1.5] :offset [4.5 0.5] :magnify 2.75)]
  (println l))

(ns cube
  (:use [geometry :only (add-vertex add-face into-faceset-factory parrallel?)])
  (:use [incanter.core :only ($= bind-columns trans)])
  (:require jreality)
  (:require planar-dual)
  (:use [planar-dual :only (M)])
  (:use [clojure.contrib.generic.math-functions :only (approx=)])
  (:use [geometry.shapes.cube])
  (:use [clojure.contrib.duck-streams :only (with-out-writer)]))

(def ifsf (if (bound? #'ifsf)
            ifsf
            (doto (de.jreality.geometry.IndexedFaceSetFactory.)
              (.setGenerateFaceNormals true))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; flags
(def *reset-geom* true)
(def *generate-nodes* true)
(def *use-backup* true)
(def *intersections* false)
(def *solution* 0)


;; legacy
(def *test-e*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Helper Fn

(defn update-geom []
  (into-faceset-factory ifsf))

(defn -main []
  (create-cube)
  (update-geom)
  (jreality/show (jreality/sgc ifsf)))

(defn reset []
  (jreality/dispose)
  (geometry/reset-geom)
  (-main))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Aux Fn




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Demo Fn
(defn demo-cc []
  (create-cube)
  (update-geom)
  (dotimes [_ 4]
    (sleep)
    (catmull-clark-subd)))

(defn demo-ds []
  (create-cube)
  (update-geom)
  (dotimes [_ 4]
    (sleep)
    (doo-sabin-subd)))

(defn demo-p []
  (create-cube)
  (update-geom)
  (dotimes [_ 2]
    (sleep)
    (planar-subd)
    (sleep)
    (sleep)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; planar stuff.

(def *test-e* '(0 0 -2/3))
(def *test-nodes-ellip*  '((-2 -2 0) (-1 -2 1) (1 -2 1)
                           (2 -2 0) (2 -1 1) (2 1 1)
                           (2 2 0) (1 2 1) (-1 2 1)
                           (-2 2 0) (-2 1 1) (-2 -1 1)))
(def *test-nodes-hyp* '((-2 -2 3) (-1 -2 2) (1 -2 1)
                        (2 -2 0) (2 -1 1) (2 1 2)
                        (2 2 3) (1 2 2) (-1 2 1)
                        (-2 2 0) (-2 1 1) (-2 -1 2)))

(def *test-nodes* *test-nodes-ellip*)

(defn show-nodes []
  (let [nds (map #(apply add-vertex %) *test-nodes*)]
    (doseq [[a b] (face-to-edges nds)]
      (geometry/add-edge a b))
    nds))

(defn show-intersection [[a b c] e]
  (let [h (geometry/plane a b c)
        dir (geometry/unit-vec ($= (geometry/skew-mat h) <*> e))
        pt  (incanter.core/solve (incanter.core/matrix [e h dir]) (incanter.core/matrix [-1 -1 0]))]
    (geometry/add-edge
     (apply add-vertex ($= pt - 2 * dir))
     (apply add-vertex ($= pt + 2 * dir)))))

(defn show-all-intersections [e pts]
  (doseq [d '(2 5 8 11)]
    (show-intersection (take 3 (drop d (cycle pts))) e)))

(defn test-intersections []
  (show-all-intersections *test-e* *test-nodes*))

(defn test-planar []
  (let [E *test-e*
        boundary *test-nodes*]
    (planar-dual/find-first E boundary)))

(defn test-qs [tp E q1]
  (let [q2 (planar-dual/project-onto-plane ($= (M (:c tp)) <*> (M (:b tp)) <*> q1) E)
        q3 (planar-dual/project-onto-plane ($= (M (:i tp)) <*> (M (:f tp)) <*> q2) E)
        q4 (planar-dual/project-onto-plane ($= (M (:g tp)) <*> (M (:h tp)) <*> q3) E)
        q5 (planar-dual/project-onto-plane ($= (M (:a tp)) <*> (M (:d tp)) <*> q4) E)]
    [q2 q3 q4 q5]))

(defn q4->q1 [tp E q4]
  (planar-dual/project-onto-plane ($= (M (:a tp)) <*> (M (:d tp)) <*> q4) E))

(defn q3->q1 [tp E q3]
  (q4->q1 tp E (planar-dual/project-onto-plane ($= (M (:g tp)) <*> (M (:h tp)) <*> q3) E)))

(defn q2->q1 [tp E q2]
  (q3->q1 tp E (planar-dual/project-onto-plane ($= (M (:i tp)) <*> (M (:f tp)) <*> q2) E)))

(declare test-qs)

(defn demo-planar-1 []
  (binding [*use-backup* false *solution* 0] (run-test))
  (sleep)
  (sleep)
  (binding [*use-backup* false *solution* 1] (run-test))
  (sleep)
  (sleep)
  (binding [*use-backup* true] (run-test)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; second planar test
(defn test-planar-2 []
  (geometry/reset-geom)
  (create-cube)
  (let [ps (planar-subd-edge-points)]
    (planar-subd-adjust-vertices ps)
    (def *test-e* (planar-subd-compute-E ps (first (geometry/get-faces))))
    (def *test-nodes* (doall (map geometry/get-vertex (planar-subd-face-boundary ps (first (geometry/get-faces))))))
    (show-nodes)))

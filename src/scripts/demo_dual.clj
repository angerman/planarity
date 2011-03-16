(ns scripts.demo-dual
  (:use [geometry :only (with-geometry create-geometry add-vertex set-vertex get-vertex add-edge into-faceset-factory
                          remove-vertex remove-node remove-edge remove-face
                          get-nodes get-edges get-faces)])
  (:use [geometry.utils :only (centroid unit-vec skew-mat edges-for-vertex faces-for-vertex)])
  (:use [incanter.core :only ($= trans)])
  (:use [geometry.subdivision.planar :only (planar-subd-compute-backup-range-pure planar-subd-compute-qs
                                                                                  select-fixed-point planar-subd-add-points-and-faces)])
  (:require [jreality :as jr])
  (:require planar-dual))

;; This code provieds a jReality window with a draggable boundary
;; Upon dragging it tries to compute a planar 3x3 patch solution
;; within the boundary.

;; The geometry structures numerates node starting form 1. While
;; the jReality backing structure starts at 0. This should be
;; kept in mind when reading the code.



(defn create-boundary []
  (let [geom (create-geometry)]
    (with-geometry geom
      (let [p1 (add-vertex -3 -3 1)
            p2 (add-vertex -1 -3 1)
            p3 (add-vertex  1 -3 1)
            p4 (add-vertex  3 -3 1)
            p5 (add-vertex  3 -1 1)
            p6 (add-vertex  3  1 1)
            p7 (add-vertex  3  3 1)
            p8 (add-vertex  1  3 1)
            p9 (add-vertex -1  3 1)
            p10 (add-vertex -3 3 1)
            p11 (add-vertex -3 1 1)
            p12 (add-vertex -3 -1 1)]
        (doseq [[a b] (partition 2 1 [p1] [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12])]
          (add-edge a b))))
    geom))

;; this code is partially adapted from the subdivision.planar stuff.

(def *corner-weight* -2/2)
(def *edge-weight* 4/2)

(defn compute-E-offset []
  "computes the offset for the central face E
   using the centroids of the corner points p1,p4,p7,p10
   and the edges p2,p3 p5,p6 p8,p9 and p11,p12

  asumes to be run within (with-geometry ...)"
  (let [c-centroid (centroid (map get-vertex [1 4 7 10]))
        e-centroid (centroid (map get-vertex [2 3 5 6 8 9 11 12]))]
    ($= *edge-weight* * e-centroid + *corner-weight* * c-centroid)))


(defn compute-E-direction []
  "the direction of E is computed as the average of the normals of the
   four faces A, C, I, G"
  (let [[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12] (map get-vertex (range 1 13))]
    (unit-vec
     (centroid (unit-vec ($= (skew-mat ($= p1 - p12)) <*> ( p2 - p1)))
               (unit-vec ($= (skew-mat ($= p4 - p3))  <*> ( p5 - p4)))
               (unit-vec ($= (skew-mat ($= p7 - p6))  <*> ( p8 - p7)))
               (unit-vec ($= (skew-mat ($= p10 - p9)) <*> ( p11 - p10)))))))

(defn compute-E []
  "Computes the supporting hyperplane for E"
  (let [offset (compute-E-offset)
        orient (compute-E-direction)]
    ($= ( -1 / ( ( trans orient) <*> offset ) ) * orient)))

(defn compute-qs []
  (let [E (compute-E)
        boundary (range 1 13)
        vboundary (map get-vertex boundary)]
    ;; the find-first fucnction has a very bad name
    ;; but it creates a structure that contains q1
    ;; as well as the line a to i, therefore it is
    ;; sufficient to compute all q1 to q4.
    ;;
    ;; it also probides an :equal flag in the meta
    ;; data of q1 signaling the quality of the
    ;; eigenvalues.
    ;;
    ;; it might fail though.
    (if-let [struct (planar-dual/find-first E vboundary)]
      ;; if the eigenvalues are equal, that implies that the soution
      ;; does not consist of two distinct fixed-points but is instead
      ;; the whole line joing both fixed-points.
      ;;
      ;; if they are not equal, we will ensure that the fixed points
      ;; are contained in the line of intersection a = A \cap E.
      (let [q1-options (if (:equal (meta (:q1 struct)))
                         (planar-subd-compute-backup-range-pure vboundary E)
                         (filter #($= (trans (:a struct)) <*> %) (:q1 struct)))]
        ;; now there are possibly multiple configurations for q1 to
        ;; q4. We will gather all of them and then use some quality
        ;; function to determine which is the ``beste'' solution.
        (let [configurations (filter (complement nil?) (map (partial planar-subd-compute-qs struct E) q1-options))]
          ;; with the set of configurations we will use the
          ;; select-fixed-point function to select the best solution.
          ;;
          ;; might fail and return nil, if non was suitable.
          (if-let [configuration (select-fixed-point vboundary configurations)]
            (planar-subd-add-points-and-faces boundary (cons (last configuration)
                                                             (butlast configuration)))
            (println "[ E R R O R ] No suitable configuration!"))))
      (println "[ E R R O R ] Could not find q1!"))))

(defn remove-solution []
  "Removes the inner 4 vertices, adjascent edges and faces and"
  (let [nodes (filter #(> % 12) (get-nodes))]
    (doseq [n nodes]
      (remove-vertex n)
      (remove-node n)
      (doseq [e (edges-for-vertex (get-edges) n)]
        (remove-edge e))
      (doseq [f (faces-for-vertex (get-faces) n)]
        (remove-face f)))))

(defn -main []
  (let [geom (create-boundary)
        isfs (jr/indexed-face-set-factory)
        update (fn [] (with-geometry geom (into-faceset-factory isfs)))
        drag (fn [_ i [a b c]]
               (if (and (<= 0 i)
                        (< i 13))
                 (do
                   (with-geometry geom
                     (set-vertex (inc i) [a b c])
                     (remove-solution)
                     (.shrink-counter geom)
                     (compute-qs))                 
                   (update)))
               (println (format "Dragged: %d." i)))]
    (update)
    (with-geometry geom
      (compute-qs)
      (update))
    (jr/show (jr/add-drag-fn (jr/sgc isfs) drag))))

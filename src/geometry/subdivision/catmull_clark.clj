(ns geometry.subdivision.catmull-clark
  (:require [geometry :as g])
  (:use [incanter.core :only ($= bind-columns)])
  (:use [geometry :only (add-vertex get-vertex add-edge add-face)])
  (:use [geometry.utils :only (centroid face-to-edges edge)]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Catmull-Clark Subdivision

(defn catmull-clark-subd-face-points []
  (loop [faces (g/get-faces)
         edges #{}
         fpts  '()
         p->fp '()
         e->fp '()
         f->fp  {}]
    (let [face  (first faces)
          ps    (map geometry/get-vertex face)
          fp    (apply add-vertex (apply centroid ps))
          p->fp (concat p->fp (partition 2 (interleave face (repeat fp))))
          edges (apply conj edges (face-to-edges face))
          e->fp (concat e->fp (partition 2 (interleave (face-to-edges face) (repeat fp))))
          fpts  (conj fpts fp)
          f->fp (assoc f->fp face fp)]
      (if (empty? (rest faces))
        {:edges edges
         :p->fp p->fp
         :e->fp e->fp
         :f->fp f->fp
         :fpts  fpts}
        (recur (rest faces) edges fpts p->fp e->fp f->fp)))))


(defn catmull-clark-subd-edge-points [struct]
  (let [fp-for-edge (fn [e] (map second (filter #(= (first %) e) (:e->fp struct))))]
    (loop [edges (:edges struct)
           edge-points '()
           p->ep '()
           e->ep  {}]
      (let [edge (first edges)
            fpts (fp-for-edge edge)
            ep (apply add-vertex (apply centroid (map geometry/get-vertex (concat (first edges) fpts))))
            edge-points (conj edge-points ep)
            p->ep (concat p->ep (partition 2 (interleave edge (repeat ep))))
            e->ep (assoc e->ep edge ep)]
        (if (empty? (rest edges))
          (merge struct
                 {:edge-points edge-points
                  :p->ep       p->ep
                  :e->ep       e->ep})
          (recur (rest edges) edge-points p->ep e->ep))))))

(defn catmull-clark-subd-vertices [struct]
  (let [pts (apply disj (set (geometry/get-nodes)) (concat (:edge-points struct)
                                                       (:fpts struct)))
        fp-for-point (fn [p] (map second (filter #(= (first %) p) (:p->fp struct))))
        ep-for-point (fn [p] (map second (filter #(= (first %) p) (:p->ep struct))))]
    (doseq [p pts]
      (let [fpts (fp-for-point p)
            epts (ep-for-point p)
            n    (count fpts)]
        (if (= n (count epts))
          (let [pt-weight (/ (- n 2) n) #_(if false (/ (- n 3) n) (- 1 (/ 7 (* 4 n))))
                ep-weight (/ 1 n) #_(if false (/ 2 n) (/ 3 (* 2 n)))
                fp-weight (/ 1 n) #_(if false (/ 1 n) (/ 1 (* 4 n)))
                pt        (geometry/get-vertex p)
                eps       (apply centroid (map geometry/get-vertex epts))
                fps       (apply centroid (map geometry/get-vertex fpts))]
            (geometry/set-vertex p 
             ($= (bind-columns pt eps fps) <*> [ pt-weight ep-weight fp-weight ])))
          (geometry/set-vertex p
           (apply centroid (geometry/get-vertex p) (map geometry/get-vertex epts))))))
    struct))

(defn catmull-clark-subd-faces [struct]
  (let [faces (geometry/get-faces)]
    (doseq [face faces]
      (doseq [[a b c] (partition 3 (interleave face
                                               (drop 1 (cycle face))
                                               (drop 2 (cycle face))))]
        (add-face b
                  (get (:e->ep struct) (edge a b))
                  (get (:f->fp struct) face)
                  (get (:e->ep struct) (edge b c))))
      (apply geometry/remove-face face))))

(defn catmull-clark-subd []
  (-> (catmull-clark-subd-face-points)
      (catmull-clark-subd-edge-points)
      (catmull-clark-subd-vertices)
      (catmull-clark-subd-faces))
  ;; remove old edges
  (doseq [edge (geometry/get-edges)]
    (apply geometry/remove-edge edge))
  ;; create new edges from faces.
  (geometry/generate-edges-from-faces))

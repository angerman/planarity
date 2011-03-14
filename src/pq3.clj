(ns pq3
  (:use [jreality])
  (:use [voss])
  (:use [planarity]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; voss surface
;;
;; The idea is the following: given u and v lines, we build a mesh of these.

;; The following specify the points in 
(def u-spec '((0.2 0.3 1)
              (0.3 0.4 1)
              (0.4 0.5 1)
              (0.5 0.4 1)
              (0.4 0.3 1)
              (0.3 0.2 1)))
(def v-spec '((0.2 0.3 1)
              (0.3 0.4 1)
              (0.4 0.5 1)
              (0.5 0.4 1)
              (0.4 0.3 1)
              (0.3 0.2 1)))

(defn knod-drag [mesh updatefn _ i pos]
  (let [nu (count @mesh)
        nv (count @mesh)
        v (mod i nv)
        u (/ (- i v) nu)]
    (if (and (> u 1) (> v 1) (< u (- nu 2)) (< v (- nv 2)))
      (do
        (dosync
         (alter mesh assoc-in [u v] pos))
        (make-primal-mesh-planar mesh u v)
        (updatefn mesh)))))

(defn update-mesh [qf mesh]
  (prn "Updading cooridnates")
  (. qf setVertexCoordinates (lst-to-double-array3 @mesh))
  (. qf update))

(defn -main []
  (prn "Here we go with Voss!")
  (let [mesh (setup-mesh [0 0 0] u-spec v-spec)]
    (do (fill-mesh mesh))
    (let [-mesh (quad-mesh @mesh)
          update (partial update-mesh -mesh)]
      (show (let [cmp (sgc nil)]
              (.addChild cmp (add-drag-fn (sgc -mesh)
                                          (partial knod-drag mesh update)))
              cmp)))))

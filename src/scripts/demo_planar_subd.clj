(ns scripts.demo-planar-subd
  (:require [geometry :as g])
  (:require [jreality :as jr])
  (:require [geometry.shapes.cube :as cube])
  (:require [geometry.shapes.doughnut :as donut])
  (:require [geometry.shapes.doughnut-2 :as donut-2])
  (:require [geometry.shapes.multi-doughnut :as multi-donut])
  (:require [geometry.shapes.star :as star])
  (:require [geometry.shapes.T :as T])
  (:require [geometry.subdivision.planar :as subd])
  (:use [utils :only (sleep)]))

(defn -main []
  (let [ifsf (jr/indexed-face-set-factory)
        update (fn [] (g/into-faceset-factory ifsf))]
    (jr/show (jr/sgc ifsf))
    (doseq [shape [cube/create donut/create donut-2/create multi-donut/create star/create T/create]]
      (g/with-geometry (shape)
        (update)
        (dotimes [_ 1]
          (sleep)
          (subd/planar-subd)
          (update))
        (sleep 5000)))))

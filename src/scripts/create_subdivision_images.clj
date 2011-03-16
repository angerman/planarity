(ns scripts.create-subdivision-images
  (:require [geometry :as g])
  (:require [geometry.renderer :as r])
  (:require [geometry.shapes.cube :as cube])
  (:require [geometry.shapes.doughnut :as doughnut])
  (:require [geometry.shapes.multi-doughnut :as mdoughnut])
  (:require [geometry.shapes.T :as T])
  (:use [geometry.subdivision.catmull-clark :only (catmull-clark-subd)])
  (:use [geometry.subdivision.doo-sabin :only (doo-sabin-subd)])
  (:use [geometry.subdivision.planar :only (planar-subd)])
  (:require [tikz :as t])
  (:require [clojure.contrib.duck-streams :as ds]))

(defn render-to-file [file geom]
  (let [f (t/resolve-file (t/temp-file file))
        eye  [5 6 4]
        up   [0 0 1]
        near 1.0
        far  13.0]
    (println (format "Writing file %s" f))
    (ds/with-out-writer f
      (println (t/header))
      (doseq [line (r/render-geometry geom eye up near far :magnify 15.0)]
        (println line))
      (println (t/footer)))))

(defn -main []
  (let [tempf (t/temp-file "temp.tex")]
    (println (format "Will use file: %s." tempf))
    (let [files
          (for [[algorithm aname] [[catmull-clark-subd 'catmull-clark-subd]
                                   [doo-sabin-subd 'doo-sabin-subd]
                                   [ planar-subd 'planar-subd]]
                [shape sname]     [[(cube/create) 'cube]
                                   [(doughnut/create) 'doughnut]
                                   [(mdoughnut/create) 'mdoughnut]
                                   [(T/create) 'T]]]
            (for [i (range 3)]
              (let [fname (str sname "." aname ".step-" i ".tikz")]
                (do (if (> i 0) (g/with-geometry shape (algorithm))))
                (render-to-file fname shape)
                (format "\\input{%s};\n" fname))))]
      (ds/with-out-writer (t/resolve-file tempf)
        (t/template (reduce str (apply concat files)))))
    (t/compile-and-show (t/resolve-file tempf))))

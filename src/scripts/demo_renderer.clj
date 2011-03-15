(ns scripts.demo-renderer
  (:require [geometry :as g])
  (:require [geometry.renderer :as r])
  (:require [geometry.shapes.cube :as cube])
  (:require [geometry.subdivision.catmull-clark :as sd])
  (:require [tikz :as t])
  (:require [clojure.contrib.duck-streams :as ds]))

(defn -main []
  (println (format "Will use file: %s." (t/temp-file "temp.tex")))
  (let [f (t/resolve-file (t/temp-file "temp.tex"))
        my-cube (cube/create)
        eye     [5 6 2.2]
        up      [0 0 1]
        near    1.0
        far     10.0]
    (let [figures (for [i (range 5)]
                    (do (if (> i 0) (g/with-geometry my-cube
                                      (sd/catmull-clark-subd)))
                        (concat [(t/header)]
                                (r/render-geometry my-cube eye up near far :magnify 15.0)
                                [(t/footer)])))]
      (ds/with-out-writer f
        (t/template
         (reduce str (apply concat figures)))))
    (t/compile-and-show f)))

(ns scripts.demo-render-schwarz-lantern
  (:require [geometry :as g])
  (:require [geometry.renderer :as r])
  (:require [geometry.shapes.schwarz-lantern :as lantern])
  (:require [tikz :as t])
  (:require [clojure.contrib.duck-streams :as ds]))

(defn -main []
  (println (format "Will use file: %s." (t/temp-file "lantern.tex")))
  (let [f (t/resolve-file (t/temp-file "lantern.tex"))
        f2 (t/resolve-file (t/temp-file "lantern2.tex"))
        l (lantern/create 5 6)
        l2 (lantern/create 10 12)
        eye [5 6 2.2]
        up  [0 0 1]
        near 1.0
        far  10.0]
    (ds/with-out-writer f
      (t/template
       (reduce str (concat (interpose "\n" (concat [(t/header)] (r/render-geometry l eye up near far :magnify 15.0) [(t/footer)]))))))
    (ds/with-out-writer f2
      (t/template
       (reduce str (concat (interpose "\n" (concat [(t/header)] (r/render-geometry l2 eye up near far :magnify 15.0) [(t/footer)]))))))
    (t/compile-and-show f)
    (t/compile-and-show f2)))

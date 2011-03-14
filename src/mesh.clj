(ns mesh
  (:use [incanter.core :only ($=)])
  (:use [geometry.utils :onlye (line-meet)]))

(defn create-mesh [nu nv]
  (ref (vec (map vec (partition nv (repeat (* nu nv)
                                           nil))))))

(def *mesh-ref*)
(declare val)
(defn get [u v] (get-in @*mesh-ref* [u v]))
(defn set [u v val] (dosync (alter *mesh-ref* assoc-in [u v] val)) val)
(defn du+ [u v] ($= (val (inc u) v) - (val u v)))
(defn du- [u v] ($= (val (dec u) v) - (val u v)))
(defn dv+ [u v] ($= (val u (inc v)) - (val u v)))
(defn dv- [u v] ($= (val u (dec v)) - (val u v)))
(defn val [u v]
  (let [n (get u v)]
    (if-not (nil? n)
      n
      (set u v (line-meet (val (dec u) v) (du- u (dec v))
                          (val u (dec v)) (dv- (dec u) v))))))

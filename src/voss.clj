(ns voss
  (:use [incanter.core :only ($=)])
  (:use [geometry])
  (:use [homogenous-geometry])
  (:use [mesh]))

(defn points-from-spec [spec]
  (loop [m id
         u spec
         acc []]
    (if-not (empty? u)
      (recur ($= m <*> (apply move (first u)))
             (rest u)
             (concat acc [($= m <*> [1 1 0 0])]))
      acc)))

(defn set-origin [mesh pt]
  (dosync
   (alter mesh assoc-in [0 0] pt)))

(defn set-u-boundary [mesh pts]
  (dosync
   (doseq [u (range 1 (count @mesh))]
     (alter mesh assoc-in [u 0] (nth pts u)))))

(defn set-v-boundary [mesh pts]
  (dosync
   (doseq [v (range 1 (count (first @mesh)))]
     (alter mesh assoc-in [0 v] (nth pts v)))))

(defn fill-mesh [mesh]
  (binding [*mesh-ref* mesh]
    (val (dec (count @mesh))
         (dec (count (first @mesh))))))

(defn setup-mesh [origin u-spec v-spec]
  (let [u-pts (map #($= (rest %) / (first %)) (points-from-spec u-spec))
        v-pts (map #($= (rest %) / (first %))
                   (map #($= (rot-z (* 2 (/ Math/PI 4))) <*> %) (points-from-spec v-spec)))
        mesh  (create-mesh (count u-pts) (count v-pts))]
    (set-origin mesh origin)
    (set-u-boundary mesh u-pts)
    (set-v-boundary mesh v-pts)
    mesh))


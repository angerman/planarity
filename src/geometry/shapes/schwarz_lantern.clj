(ns geometry.shapes.schwarz-lantern
  (:use [geometry :only (create-geometry with-geometry add-vertex add-face generate-edges-from-faces)]))

(defn create-ring [m z odd]
  "Creates a set of m evenly spaced vertices lying on a circle at level z.
   odd controls if the vertices are adjusted by PI/(2m)."
  (for [angle (map #(/ (+ % (if odd 1/2 0)) m) (range m))]
    (add-vertex (Math/cos (* 2 Math/PI angle))
                (Math/sin (* 2 Math/PI angle))
                z)))

(defn create-faces [ring-a ring-b]
  "Creates trinagle faces between both rings."
  (doseq [face (partition 3 1 [(first ring-a)] (interleave ring-a ring-b))]
    (apply add-face face))
  (add-face (last ring-b) (first ring-a) (first ring-b)))

(defn create [m n]
  "Create geometry for a Schwarz Lantern of n rings
   with each m vertices."
  (let [lantern (create-geometry)]
    (with-geometry lantern
      (let [r (for [i (range n)] (create-ring m (* i (/ 3 n)) (odd? i)))]
        (loop [rings (rest r)
               ring  (first r)
               swap  false]
          (if-not (empty? rings)
            (do
              (if swap
                (create-faces (first rings) ring)
                (create-faces ring (first rings)))
              (recur (rest rings)
                     (first rings)
                     (not swap))))))
      (generate-edges-from-faces))
    lantern))

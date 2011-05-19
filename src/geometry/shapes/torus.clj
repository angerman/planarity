(ns geometry.shapes.torus
  (:use [geometry :only (create-geometry with-geometry add-vertex add-face generate-edges-from-faces)]))

(defn create-rings [m n r R]
  ""
  (doall
   (for [angle  (map #(/ % m) (range m))
         angle2 (map #(/ % n) (range n))]
     (add-vertex (* (Math/cos (* 2 Math/PI angle)) (+ R (* (Math/cos (* 2 Math/PI angle2)) r)))
                 (* (Math/sin (* 2 Math/PI angle)) (+ R (* (Math/cos (* 2 Math/PI angle2)) r)))
                 (* (Math/sin (* 2 Math/PI angle2)) r)))))

(defn create [m n r R]
  "Create geometry for a Torus of n rings
   with each m vertices."
  (let [torus (create-geometry)]
    (with-geometry torus
      (let [V (create-rings m n r R)
            rings (partition n V)]
        (doall
         (for [[A B] (partition 2 1 [(first rings)] rings)]
           (doall
            (for [[[a b] [c d]]
                  (let [vs (partition 2 (interleave A B))]
                    (partition 2 1 [(first vs)] vs))]
              (add-face a b d c))))))
      (generate-edges-from-faces))
    torus))

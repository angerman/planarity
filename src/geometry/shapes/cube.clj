(ns geometry.shapes.cube
  (:use [geometry :only (create-geometry with-geometry add-vertex add-face generate-edges-from-faces)] ))

(defn create []
  (let [shape (create-geometry)]
    (with-geometry shape
      (let [a (add-vertex 0 0 0)
            b (add-vertex 3 0 0)
            c (add-vertex 3 3 0)
            d (add-vertex 0 3 0)
            e (add-vertex 0 0 3)
            f (add-vertex 3 0 3)
            g (add-vertex 3 3 3)
            h (add-vertex 0 3 3)]
        (add-face a b c d)
        (add-face h g f e)
        (add-face a e f b)
        (add-face b f g c)
        (add-face c g h d)
        (add-face d h e a)
        (generate-edges-from-faces)))
    shape))

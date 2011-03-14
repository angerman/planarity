(ns geometry.shapes.T
  (:use [geometry :only (create-geometry with-geometry add-vertex add-face generate-edges-from-faces)]))

(defn create []
  (let [shape (create-geometry)]
    (with-geometry shape
      (let [a (add-vertex 1 1 0)
            b (add-vertex 1 2 0)
            c (add-vertex 1 2 1)
            d (add-vertex 1 3 1)
            e (add-vertex 1 3 2)
            f (add-vertex 1 2 2)
            g (add-vertex 1 1 2)
            h (add-vertex 1 0 2)
            i (add-vertex 1 0 1)
            j (add-vertex 1 1 1)

            a* (add-vertex 0 1 0)
            b* (add-vertex 0 2 0)
            c* (add-vertex 0 2 1)
            d* (add-vertex 0 3 1)
            e* (add-vertex 0 3 2)
            f* (add-vertex 0 2 2)
            g* (add-vertex 0 1 2)
            h* (add-vertex 0 0 2)
            i* (add-vertex 0 0 1)
            j* (add-vertex 0 1 1)]
        (add-face a b c j)
        (add-face i j g h)
        (add-face j c f g)
        (add-face c d e f)

        (add-face a* j* c* b*)
        (add-face i* h* g* j*)
        (add-face j* g* f* c*)
        (add-face c* f* e* d*)

        (add-face a b b* a*)
        (add-face b b* c* c)
        (add-face c c* d* d)
        (add-face d d* e* e)
        (add-face e e* f* f)
        (add-face f f* g* g)
        (add-face g g* h* h)
        (add-face h h* i* i)
        (add-face i i* j* j)
        (add-face j j* a* a)
        (generate-edges-from-faces)))
    shape))

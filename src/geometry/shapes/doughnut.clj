(ns geometry.shapes.doughnut
  (:use [geometry :only (create-geometry with-geometry add-vertex add-face generate-edges-from-faces)]))

(defn create []
  (let [shape (create-geometry)]
    (with-geometry shape
      (let [a (add-vertex 0 0 0)
            b (add-vertex 3 0 0)
            c (add-vertex 3 3 0)
            d (add-vertex 0 3 0)
            e (add-vertex 1 1 0)
            f (add-vertex 2 1 0)
            g (add-vertex 2 2 0)
            h (add-vertex 1 2 0)

            a* (add-vertex 0 0 3)
            b* (add-vertex 3 0 3)
            c* (add-vertex 3 3 3)
            d* (add-vertex 0 3 3)
            e* (add-vertex 1 1 3)
            f* (add-vertex 2 1 3)
            g* (add-vertex 2 2 3)
            h* (add-vertex 1 2 3)]
        ;; bottom
        (add-face a e f b)
        (add-face b f g c)
        (add-face c g h d)
        (add-face d h e a)
        ;; top
        (add-face a* b* f* e*)
        (add-face b* c* g* f*)
        (add-face c* d* h* g*)
        (add-face d* a* e* h*)
        ;; inner shell
        (add-face f e e* f*)
        (add-face g f f* g*)
        (add-face h g g* h*)
        (add-face e h h* e*)
        ;; outer shell
        (add-face a b b* a*)
        (add-face b c c* b*)
        (add-face c d d* c*)
        (add-face d a a* d*)
        (generate-edges-from-faces)))
    shape))

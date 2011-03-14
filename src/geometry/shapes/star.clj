(ns geometry.shapes.star
  (:use [geometry :only (create-geometry with-geometry add-vertex add-face generate-edges-from-faces)]))

(defn create []
  (let [shape (create-geometry)]
    (with-geometry shape
      (let [m (add-vertex 1 1 0)
            n (add-vertex 2 1 0)
            o (add-vertex 2 2 0)
            p (add-vertex 1 2 0)

            b- (add-vertex 1 0 1)
            c- (add-vertex 2 0 1)
            e- (add-vertex 3 1 1)
            f- (add-vertex 3 2 1)
            h- (add-vertex 2 3 1)
            i- (add-vertex 1 3 1)
            k- (add-vertex 0 2 1)
            l- (add-vertex 0 1 1)
            
            m- (add-vertex 1 1 1)
            n- (add-vertex 2 1 1)
            o- (add-vertex 2 2 1)
            p- (add-vertex 1 2 1)

            b+ (add-vertex 1 0 2)
            c+ (add-vertex 2 0 2)
            e+ (add-vertex 3 1 2)
            f+ (add-vertex 3 2 2)
            h+ (add-vertex 2 3 2)
            i+ (add-vertex 1 3 2)
            k+ (add-vertex 0 2 2)
            l+ (add-vertex 0 1 2)
            
            m+ (add-vertex 1 1 2)
            n+ (add-vertex 2 1 2)
            o+ (add-vertex 2 2 2)
            p+ (add-vertex 1 2 2)

            m* (add-vertex 1 1 3)
            n* (add-vertex 2 1 3)
            o* (add-vertex 2 2 3)
            p* (add-vertex 1 2 3)]
        ;; bottom
        (add-face n m p o)
        ;; first ring of walls    
        (add-face n n- m- m)
        (add-face o o- n- n)
        (add-face p p- o- o)
        (add-face m m- p- p)
        ;; first floor
        (add-face n- c- b- m-)
        (add-face o- f- e- n-)
        (add-face p- i- h- o-)
        (add-face m- l- k- p-)
        ;; second ring of walls
        (add-face c- c+ b+ b-)
        (add-face n- n+ c+ c-)
        (add-face e- e+ n+ n-)
        (add-face f- f+ e+ e-)
        (add-face o- o+ f+ f-)
        (add-face h- h+ o+ o-)
        (add-face i- i+ h+ h-)
        (add-face p- p+ i+ i-)
        (add-face k- k+ p+ p-)
        (add-face l- l+ k+ k-)
        (add-face m- m+ l+ l-)
        (add-face b- b+ m+ m-)
        ;; second floor
        (add-face n+ m+ b+ c+)
        (add-face o+ n+ e+ f+)
        (add-face p+ o+ h+ i+)
        (add-face m+ p+ k+ l+)
        ;; third ring of walls
        (add-face n+ n* m* m+)
        (add-face o+ o* n* n+)
        (add-face p+ p* o* o+)
        (add-face m+ m* p* p+)
        ;;top
        (add-face n* o* p* m*)
        (generate-edges-from-faces)))
    shape))

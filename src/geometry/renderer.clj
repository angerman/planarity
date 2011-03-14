(ns geometry.render
  (:use [incanter.core :only ($= matrix)]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Render functions (-> TikZ)

(defn render-node-z-weight [node]
  (/ (geometry/norm (geometry/get-vertex node))
     (Math/sqrt 3)))

(defn render-edge-weight [edge]
  (* 1/2 (reduce + (map render-node-z-weight edge))))

(defn render-face-weight [face]
  (/ (reduce + (map render-node-z-weight face)) (count face)))

;; target range 25 -- 100
(defn render-depth-value [depth max-depth]
  (+ 25 (* 75 (- 1 (/ depth max-depth)))))


;; comparators
(defn render-comp-face [face-a face-b]
  (> (render-face-weight face-a)
     (render-face-weight face-b)))

(defn render-comp-edge [edge-a edge-b]
  (> (render-edge-weight edge-a)
     (render-edge-weight edge-b)))

(defn render-vertex-to-2d [vertex]
  (let [offset [1.5 1.5]
        scale  1
        m (matrix [[0.867 -0.9396  0.0]
                   [0.5    0.3402 -1.0]])]
    ($= offset + scale * m <*> vertex)))

(defn render-nodes-to-tikz-coordinates [nodes]
  (doseq [node nodes]
    (let [vertex (geometry/get-vertex node)
          [x y] (render-vertex-to-2d vertex)]
      (print (format "\\coordinate (node-%d) at (%f,%f);\n" node x y)))))

(defn render-nodes-to-tikz-fill-circles [nodes]
  (let [snodes (sort #(< (render-node-z-weight %) (render-node-z-weight %2)) nodes)
        n->depth (apply hash-map (interleave snodes (map render-node-z-weight snodes)))
        max-depth (apply max (vals n->depth))]
    (print "\\foreach \\n/\\w in {")
    (print
     (reduce #(str % ", " %2)
             (for [node nodes]
               (let [weight (render-depth-value (get n->depth node) max-depth)]
                 (format "node-%d/%3.5f" node weight)))))
    (print "}\n")
    (print (format "\t\\fill[black!\\w] (\\n) circle (2pt);\n"))))

(defn render-edges-to-tikz [edges]
  (let [sedges (sort render-comp-edge edges)
        e->depth (apply hash-map (interleave sedges (map render-edge-weight sedges)))
        max-depth (apply max (map render-node-z-weight (set (reduce concat sedges))))]
    (for [edge sedges]
      (let [weight (render-depth-value (get e->depth edge) max-depth)
            [a b] edge]
        `(~weight ~(format "\\draw[black!%3.5f] (node-%d) -- (node-%d);\n" weight a b))))))

(defn render-faces-to-tikz [faces]
  (let [sfaces (sort render-comp-face faces)
        f->depth (apply hash-map (interleave sfaces (map render-face-weight sfaces)))
        max-depth (apply max (map render-node-z-weight (set (reduce concat sfaces))))]
    (for [face sfaces]
      (let [weight (render-depth-value (get f->depth face) max-depth)]
        `(~weight ~(format "\\fill[black!%3.5f] %s -- cycle;\n" weight (reduce #(str % " -- " %2) (map #(format "(node-%d)" %) face))))))))

(defn render [& {:keys [faces] :or {faces false}}]
  (render-nodes-to-tikz-coordinates (geometry/get-nodes))
  (doseq [[_ t] (sort-by first (concat
                                (if faces (render-faces-to-tikz (geometry/get-faces)) '())
                                (render-edges-to-tikz (geometry/get-edges))))]
    (print t))
  ;; they would need to be drawn interleaved with the edges.
  ;;(render-nodes-to-tikz-fill-circles (geometry/get-nodes))
  )

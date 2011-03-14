(ns homogenous-geometry
  (:use [incanter.core :only ($= matrix cos sin)]))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 3d homogenous coordinates (+ helper functions)

(def id [[1 0 0 0]
         [0 1 0 0]
         [0 0 1 0]
         [0 0 0 1]])

(def origin [1 0 0 0])

(defn rot-x [t]
      (matrix [[1 0 0 0] [0 1 0 0] [0 0 (cos t) (* -1 (sin t))] [0 0 (sin t) (cos t)]]))
(defn rot-z [a]
      (matrix [[1 0 0 0] [0 (cos a) (* -1 (sin a)) 0] [0 (sin a) (cos a) 0] [0 0 0 1]]))
(defn translate [[x y z]]
      (matrix [[1 0 0 0] [x 1 0 0] [y 0 1 0] [z 0 0 1]]))
(defn move [t a l]
  ($= (rot-x t) <*> (rot-z a) <*> (translate [l 0 0])))


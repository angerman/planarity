;; -*- coding: utf-8 -*-
;; computing planar quadlieral nets.
;; 
;; using: JReality + Incanter + Clojure
;; 
;; Author: Moritz Angermann <moritz.angermann@gmail.com>
;; Date: May 24th 2010
;; TODO: reflect changes to dual boundary points

(import 'java.awt.Color)
(import '(java.awt.event KeyListener KeyEvent))
(import 'de.jreality.shader.ShaderUtility)
(import '(de.jreality.scene Appearance Geometry SceneGraphComponent))
(import '(de.jreality.scene.data Attribute StorageModel IntArrayArray$Array))
(import 'de.jreality.tools.DragEventTool)
(import 'de.jreality.plugin.JRViewer)
(import 'de.jreality.geometry.IndexedFaceSetFactory)
(use 'incanter.core) ;; For $=, decompose-eig and general matrix computation

;; twister matrix for the cross product
(def T1 (matrix [[0 1 0]
                 [0 0 1]
                 [1 0 0]]))
(def T2 (matrix [[0 0 1]
                 [1 0 0]
                 [0 1 0]]))
(defn cross-product [x y]
  ($= (T1 <*> x) * (T2 <*> y) - (T2 <*> x) * (T1 <*> y)))

(defn norm [x]
  (sqrt ($= (trans x) <*> x)))

;; the normal of three points.
(defn normal [a b c]
  (let [u ($= a - b)
        v ($= c - b)
        n (cross-product u v)]
    ($= n / (norm n))))

;; turn a point (3d) into homogenous coordinates
(defn ->hc
  ([[x y z]] [x y z 1])
  ([[x y z] d] [x y z d]))

;; turn a homogenous coordiante into a point in 3d space.
(defn <-hc [[x y z d]]
  ($= [x y z] / d))

;; Hyperplane in Homogenous Coordinates
(defn HCHyperplane [a b c]
  (let [n (normal a b c)
        d ($= -1 * (trans a) <*> n)]
    (->hc n d)))


;; Projection Matrix. To project any point onto a hyperplane in the direction of a given origin.
(def I4 (identity-matrix 4))
(defn M [origin hyperplane]
  ($= origin <*> (trans hyperplane) - (trans hyperplane) <*> origin * I4))


;; all the heavy lifting to find planar quadriliterals.
;; 1 ----- 2 - m - 3 ----- 4  a = X ^ m
;; |       |       |       |  b = X ^ n
;; |   I   |       |   L   |  c = X ^ o
;; |       |       |       |  d = X ^ p
;; 5 ----- 6 ----- 7 ----- 8
;; |       |p1   p4|       |  
;; n       |   X   |       p
;; |       |p2   p3|       |
;; 9 ---- 10 ---- 11 ---- 12
;; |       |       |       |
;; |   J   |       |   K   |
;; |       |       |       |
;;13 ---- 14 - o - 15 --- 16
(defn compute-points [points _ P]
  (let [pts (apply vector (map matrix points))
        x (<-hc (matrix P))
        ; xn ($= x - (matrix [0 0 1])) ;; We get P in Homogenous Coordinates
        xn ($= x / (norm x)) ;; scale x so it's always on the unit sphere
        ;; x ($= xn + (matrix [0 0 1]))
        pt (fn [n] (nth pts n)) ;; return the point at index n from pts (points)
        hcpt (comp ->hc pt) ;; return the point at index n in homogenous coordinates
        ;; Hyperplanes
        I (apply HCHyperplane (map pt '(4 0 1)))
        J (apply HCHyperplane (map pt '(13 12 8)))
        K (apply HCHyperplane (map pt '(11 15 14)))
        L (apply HCHyperplane (map pt '(2 3 7)))
        X (matrix P)
        ;; Intersection points
        a ($= (M (hcpt 1) X)  <*> (hcpt 2)) ;; intersetcion of a with X
        b ($= (M (hcpt 8) X)  <*> (hcpt 4)) ;; intersection of b with X
        c ($= (M (hcpt 14) X) <*> (hcpt 13))
        d ($= (M (hcpt 7) X)  <*> (hcpt 11))
        ;; Hyperplane projections
        Mx (M (->hc x) X)
        Mj (M b J) ;; p_2 = Mj p_1
        Mk (M c K) ;; p_3 = Mk p_2
        Ml (M d L) ;; p_4 = Ml p_3
        Mi (M a I) ;; p_1 = Mj p_4
        ;; Combined Projection. (Any point -Mx-> p1 -Mj-> p2 -Mk-> p3 -Ml-> p4 -Mi-> p1)
        ;; Eigenvectors of this are candidates for p1 (They lie in I and X)
        M+ ($= Mi <*> Ml <*> Mk <*> Mj <*> Mx)
        eig (decomp-eigenvalue M+)
        val (map abs (:values eig))
        ;; this was _really_ nasty, because iterating over the :vectors matrix returned _rows_
        ;; and not the expected vectors. Iterating over the transposed matrix and transposing
        ;; the result gives the expected eigen-vectors. A more geometric approach would be
        ;; the mapping of the unit vectors e_1, ..., e_4 over the (:vectors eig) matrix. Though
        ;; this would incure uneccessary computation.
        vs+ (map trans (trans (:vectors eig)))
        ;; find the Eigenvector corresponding to the largest Eigenvalue
        v   (second (apply max-key first (map vector val vs+)))]
    (if (not= v (matrix [0 0 0 1]))
      (let [[p1 p2 p3 p4] (reductions #($= %2 <*> %) v [Mj Mk Ml])]
        ;; create some statisitc for each eigenvector. To determin it's performance.
        ;; Will show which Eigenvector has been choosen `<-' and what it's Eigen-value is
        ;; as well as the if the point lies within X I J K or L. Points are
        ;; displayed as * (active) or _ (inactive).
        (println (format "%-5s%-5s%-5s%-5s%-5s    %s" "X" "I" "J" "K" "L" "Eigen-value"))
        (doseq [[v- v+] (map vector vs+ (:values eig))]
          (doseq [H (map trans [X I J K L])]
            (doseq [p (reductions #($= %2 <*> %) v- [Mj Mk Ml])]
              (print (if ($= (abs ($= H <*> p)) < 1E-13) "*" "_")))
            (print " "))
          (print (if ($= v == v-) "<- " "   "))
          (println (format "%+.3e" v+)))
        ;; given all four points p1,...,p4 these represent hyperplanes in their dual.
        (let [i I
              b (HCHyperplane (pt 8) (pt 4) (<-hc p1))
              j J
              c (HCHyperplane (pt 14) (pt 13) (<-hc p2))
              k K
              d (HCHyperplane (pt 7) (pt 11) (<-hc p3))
              l L
              a (HCHyperplane (pt 1) (pt 2) (<-hc p4))
              x X]
          (do
            ;; update the dragged point
            (doseq [pos '(0 1 2)]
              (aset-double points 16 pos (nth x pos)))
            ;; update primal points
            (if (empty? (filter #($= (nth % 3) == 0) [p1 p2 p3 p4]))
              (doseq [[p idx] (map vector [p1 p4 p2 p3] [5 6 9 10])
                      pos '(0 1 2)]
                (aset-double points idx pos (nth (<-hc p) pos))))
            ;; update dual points
            (doseq [[p idx] (map vector [i a l b x d j c k] (map #(+ 17 %) (range 9)))
                    pos '(0 1 2)]
              (if ($= (nth p 3) != 0)
                (aset-double points idx pos (nth (<-hc p) pos)))))))))
  points)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Scene Definition
(defn to-type-array [type coll]
  (into-array (map #(into-array type %) coll)))
(def to-double-array (partial to-type-array Double/TYPE))
(def to-int-array (partial to-type-array Integer/TYPE))

;; Vertices
(def primal-v
     [[-2,-1.5,-2],[-2,-1,-1],[-2,+1,-1],[-2,+1.5,-2],
      [-1,-2,-1],[-1,-1,+1],[-1,+1,+1],[-1,+2,-1],
      [+1,-2,-1],[+1,-1,+1],[+1,+1,+1],[+1,+2,-1],
      [+2,-1,-2],[+2,-1,-1],[+2,+1,-1],[+2,+1,-2],
      [ 0, 0, 2]])
(def dual-v
     [[0 0 0] [0 0 0] [0 0 0]
      [0 0 0] [0 0 0] [0 0 0]
      [0 0 0] [0 0 0] [0 0 0]])
(def v (to-double-array (concat primal-v dual-v)))

;; Draggable Vertices
;; generic swapping function generator
(defn swapper [a b] #(if (= a %) b a))

(def primal-movable #{0,1,2,3,4,7,8,11,12,13,14,15})
(def dual-movable #{17,19,21,23,25})
(def movable (atom primal-movable))
(defn swap-movable [] (swap! movable (swapper primal-movable dual-movable)))

(def primal-edges
     (to-int-array
      [        [ 0  1]         [ 1  2]         [ 2  3]
       [ 0  4]         [ 1  5]         [ 2  6]         [ 3  7]
               [ 4  5]         [ 5  6]         [ 6  7]
       [ 4  8]         [ 5  9]         [ 6 10]         [ 7 11]
               [ 8  9]         [ 9 10]         [10 11]
       [ 8 12]         [ 9 13]         [10 14]         [11 15]
               [12 13]         [13 14]         [14 15]
      ]))
(def dual-edges
     (to-int-array
      [        [17 18]         [18 19]
       [17 20]         [18 21]         [19 22]
               [20 21]         [21 22]
       [20 23]         [21 24]         [22 25]
               [23 24]         [24 25]
      ]))
(def e (atom primal-edges))
(defn swap-edges [] (swap! e (swapper primal-edges dual-edges)))

;; Faces
(def primal-faces
     (to-int-array
      [[ 0, 4, 5, 1],[ 1, 5, 6, 2],[ 2, 6, 7, 3],
       [ 4, 8, 9, 5],[ 5, 9,10, 6],[ 6,10,11, 7],
       [ 8,12,13, 9],[ 9,13,14,10],[10,14,15,11]]))
(def dual-faces
     (to-int-array
      [[17 20 21 18] [18 21 22 19]
       [20 23 24 21] [21 24 25 22]]))
(def f (atom primal-faces))
(defn swap-faces [] (swap! f (swapper primal-faces dual-faces)))
(def fc (into-array Color (concat [] ;(repeat 9 Color/green)
                                  (repeat 9 Color/blue))))
;; FaceSet Geometry
(def fs (IndexedFaceSetFactory.))
(doto fs
  (.setVertexCount (alength v))
  (.setVertexCoordinates v)
  (.setFaceCount (alength @f))
  (.setFaceIndices @f)
  (.setFaceColors fc)
  (.setEdgeCount (alength @e))
  (.setEdgeIndices @e)
  (.setGenerateFaceNormals true)
  (.update))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; DragTool
(defmacro update-points [f ev]
  `(let [ps#  (. ~ev getPointSet)
         n#   (. ps# getNumPoints)
         pts# (make-array Double/TYPE n# 3)
         COORDS# Attribute/COORDINATES
         STORAGEMODEL# (. StorageModel/DOUBLE_ARRAY array 3)]
     (. (. ps# getVertexAttributes COORDS#) toDoubleArrayArray pts#)
     (. ps# setVertexAttributes COORDS# (. STORAGEMODEL# createReadOnly (~f pts# (. ~ev getIndex) (. ~ev getPosition))))))

(defn just-move [points index point]
  ;; update point
  (aset points index point)
  ;; return point-set
  points)

(def dt (DragEventTool.))
(.addPointDragListener dt (proxy [de.jreality.tools.PointDragListener] []
                            (pointDragStart [e])    ;; nothing
                            (pointDragEnd [e])      ;; nothing
                            (pointDragged [e] ;; magic happens here :D
                                          (cond
                                           (contains? @movable (. e getIndex)) (update-points just-move e)
                                           (= 16 (. e getIndex))              (update-points compute-points e)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Appearance setup
(defn setup-appearance [ap]
  (let [gs (ShaderUtility/createDefaultGeometryShader ap true)
        ps (. gs createPolygonShader "default")
        ls (. gs createLineShader "default")
        pts (. gs createPointShader "default")]
    (doto gs
      (.setShowFaces true)
      (.setShowLines true)
      (.setShowPoints true))
    (doto ps
      (.setDiffuseColor Color/blue))
    (doto ls
      (.setDiffuseColor Color/yellow)
      (.setTubeRadius 0.03))
    (doto pts
      (.setDiffuseColor Color/red)
      (.setPointRadius 0.05)))
  ap)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Component Graph assembly and Viewer initialization
(def cmp (SceneGraphComponent.))
(doto cmp
  (.setGeometry (.getGeometry fs))
  (.setAppearance (setup-appearance (Appearance.)))
  (.addTool dt))

(def viewer (JRViewer/display cmp))
;; adding a listener for the "D" key, which will dualize the image. (hide faces and edges)
(.. viewer (getViewingComponent)
    (addKeyListener
     (proxy [KeyListener] []
       (keyPressed [e]) ;; ignore
       (keyTyped [e])   ;; ignore
       (keyReleased [e]
                    (if (= (. e getKeyCode) KeyEvent/VK_D)
                      (do
                        (swap-movable)
                        (doto (. cmp getGeometry)
                          (.setFaceCountAndAttributes Attribute/INDICES (IntArrayArray$Array.  (swap-faces)))
                          (.setEdgeCountAndAttributes Attribute/INDICES (IntArrayArray$Array.  (swap-edges))))
                        (println "Dualizing")))))))
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

;; (def dt (DragEventTool.))
;; (.addPointDragListener dt (proxy [de.jreality.tools.PointDragListener] []
;;                             (pointDragStart [e])    ;; nothing
;;                             (pointDragEnd [e])      ;; nothing
;;                             (pointDragged [e] ;; magic happens here :D
;;                                           (cond
;;                                            (contains? @movable (. e getIndex)) (update-points just-move e)
;;                                            (= 16 (. e getIndex))              (update-points compute-points e)))))

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
;;  (.addTool dt)
)

(def viewer (JRViewer/display cmp))
;; adding a listener for the "D" key, which will dualize the image. (hide faces and edges)
;; (.. viewer (getViewingComponent)
;;     (addKeyListener
;;      (proxy [KeyListener] []
;;        (keyPressed [e]) ;; ignore
;;        (keyTyped [e])   ;; ignore
;;        (keyReleased [e]
;;                     (if (= (. e getKeyCode) KeyEvent/VK_D)
;;                       (do
;;                         (swap-movable)
;;                         (doto (. cmp getGeometry)
;;                           (.setFaceCountAndAttributes Attribute/INDICES (IntArrayArray$Array.  (swap-faces)))
;;                           (.setEdgeCountAndAttributes Attribute/INDICES (IntArrayArray$Array.  (swap-edges))))
;;                         (println "Dualizing")))))))
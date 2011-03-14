(ns jreality
  (:import [de.jreality.plugin JRViewer])
  (:import [de.jreality.scene Appearance Geometry SceneGraphComponent])
  (:import [java.awt Color])
  (:import [java.awt.event KeyListener KeyEvent])
  (:import [de.jreality.shader ShaderUtility])
  (:import [de.jreality.geometry QuadMeshFactory IndexedLineSetFactory IndexedFaceSetFactory])
  (:import [de.jreality.tools DragEventTool]))


;; creates a default appearance
(defn default-appearance []
  (let [ap (Appearance.)
        gs (ShaderUtility/createDefaultGeometryShader ap true)
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
      (.setTubeRadius 0.02)) ;; 0.03
    (doto pts
      (.setDiffuseColor Color/red)
      (.setPointRadius 0.05)) ;; 0.08
    ap))

(defn red-trans-appearance []
  (let [ap (Appearance.)
        gs (ShaderUtility/createDefaultGeometryShader ap true)
        ps (. gs createPolygonShader "default")]
    (doto gs
      (.setShowFaces true))
    (doto ps
      (.setDiffuseColor Color/red)
      (.setTransparency 0.8))
    ap))

(defn green-trans-appearance []
  (let [ap (Appearance.)
        gs (ShaderUtility/createDefaultGeometryShader ap true)
        ps (. gs createPolygonShader "default")]
    (doto gs
      (.setShowFaces true))
    (doto ps
      (.setDiffuseColor Color/green)
      (.setTransparency 0.8))
    ap))

(defn lst-to-double-array3 [lst]
  (let [#^"[[[D" arr (make-array Double/TYPE (count lst) (count (first lst)) (count (ffirst lst)))]
    (dotimes [u (count lst)]
      (dotimes [v (count (nth lst u))]
        (aset arr (int u) (int v) #^doubles (double-array (nth (nth lst u) v)))))
    arr))

(defn lst-to-double-array2 [lst]
  (let [#^"[[D" arr (make-array Double/TYPE (count lst) (count (first lst)))]
    (dotimes [u (count lst)]
      (aset arr (int u) #^doubles (double-array (nth lst u))))
    arr))


(defn indexed-face-set-factory []
  (doto (IndexedFaceSetFactory.)
    (.setGenerateFaceNormals true)))

(defn quad-mesh [lst]
  (let [f (QuadMeshFactory.)]
    (doto f
      (.setULineCount (count lst))
      (.setVLineCount (count (first lst)))
      (.setClosedInUDirection false)
      (.setClosedInVDirection false)
      (.setVertexCoordinates (lst-to-double-array3 lst))
      (.setGenerateFaceNormals true)
      (.setGenerateTextureCoordinates true)
      (.setGenerateEdgesFromFaces true)
      (.setEdgeFromQuadMesh true)
      (.update))
    f))

(defn line-set [lst edges]
  (let [f (IndexedLineSetFactory.)]
    (doto f
      (.setVertexCount (count lst))
      (.setVertexCoordinates (lst-to-double-array2 lst))
      (.setEdgeCount (/ (count edges) 2))
      (.setEdgeIndices (int-array edges) 2)
      (.update))
    f))

(defn sgc [f & ap]
  (let [cmp (SceneGraphComponent.)]
    (if f
      (.setGeometry cmp (.getGeometry f)))
    (.setAppearance cmp (if-not (nil? ap) (first ap) (default-appearance)))
    cmp))

(def *viewer* (ref #{}))

(defn last-viewer []
  (JRViewer/getLastJRViewer))

(defn show [cmp]
  (JRViewer/display cmp)
  (dosync
   (alter *viewer* conj (last-viewer))))

(defn dispose
  ([viewer] (try
              (.dispose viewer)
              (catch Exception ex)
              (finally
               (dosync (alter *viewer* disj viewer)))))
  ([] (dispose (first @*viewer*))))

(defn drag-tool [f]
  (let [dt (DragEventTool.)]
    (.addPointDragListener dt (proxy [de.jreality.tools.PointDragListener] []
                                (pointDragStart [e])
                                (pointDragEnd [e])
                                (pointDragged [e]
                                              (let [i (.getIndex e)
                                                    x (.getX e)
                                                    y (.getY e)
                                                    z (.getZ e)]
                                                (f e i [x y z])))))
    dt))

(defn add-tool [cmp tool]
  (.addTool cmp tool)
  cmp)

(defn add-drag-fn [cmp f]
  (add-tool cmp (drag-tool f)))

(def *key-map*
  {:V KeyEvent/VK_V
   :R KeyEvent/VK_R
   :S KeyEvent/VK_S
   :M KeyEvent/VK_M})

(defn add-key-event [viewer key f]
  (.. viewer
      (getViewer)
      (getViewingComponent)
      (addKeyListener
       (proxy [KeyListener] []
         (keyPressed [e])
         (keyTyped [e])
         (keyReleased [e]
                      (if (= (. e getKeyCode) (key *key-map*))
                        (f)))))))

(ns tikz
  (:require [clojure.contrib.shell-out :as shell])
  (:require [clojure.contrib.string :as str])
  (:import (java.io File)))

(def *Skim:displayline* "/Applications/Skim.app/Contents/SharedSupport/displayline")
(def *XeLaTeX:cmd*      "/usr/texbin/xelatex")
(def *temp:dir*         "~/temp")

(defn header []
  "\\begin{tikzpicture}")

(defn footer []
  "\\end{tikzpicture}")

(defn template-header []
  ((comp (partial str/join "\n") list)
   "\\documentclass{article}"
   "\\usepackage{tikz,nicefrac,amsmath,pifont}"
   "\\usetikzlibrary{arrows,decorations,backgrounds,patterns,matrix,shapes,fit,calc,shadows,plotmarks,intersections}"
   "\\usepackage[graphics,tightpage,active]{preview}"
   "\\PreviewEnvironment{tikzpicture}"
   "\\newlength{\\imagewidth}"
   "\\newlength{\\imagescale}"
   "\\begin{document}"))

(defn template-footer []
  "\\end{document}")

(defn temp-file [file-name]
  (str/join File/separator [*temp:dir* file-name]))

(defn resolve-file [file-str]
  (let [#^String s (str file-str)
        s (.replaceAll (re-matcher #"[/\\]" s) File/separator)
        s (if (.startsWith s "~")
            (str (System/getProperty "user.home")
                 (if-not (.startsWith (subs s 1) File/separator) File/separator) (subs s 1))
            s)]
    s))

(defn get-directory [file-str]
  (->> (resolve-file file-str)
       (str/reverse)
       (str/split (re-pattern File/separator) 2)
       (last)
       (str/reverse)))

(defn template [contents]
  (doseq [line [(template-header) contents (template-footer)]]
    (println line)))

(defn compile [file-str]
  (shell/with-sh-dir (get-directory file-str)
    (shell/sh *XeLaTeX:cmd* "--interaction=nonstopmode" (resolve-file file-str) :return-map true)))

(defn show [file-str]
    (shell/sh *Skim:displayline* "-r" "0" (resolve-file file-str)))

(defn compile-and-show [file-str]
  (let [pdf-str (str/replace-re #".tex$" ".pdf" file-str)
        result  (compile file-str)]
    (if (= 0 (:exit result))
      (show pdf-str)
      (print (:out result)))))

;; Export docs/orgmode to docs/source RST.
;; ox-rst is the default. Pandoc runs only when CPMDC_ORG_EXPORT=pandoc.

(defun cpmdc-publish-with-pandoc ()
  "Export Org sources with pandoc and docs/org-links.lua."
  (let ((base (expand-file-name "orgmode" default-directory))
        (out-root (expand-file-name "source" default-directory)))
    (dolist (org-file (directory-files-recursively base "\\.org\\'"))
      (let* ((rel (file-relative-name org-file base))
             (rst-file
              (expand-file-name
               (concat (file-name-sans-extension rel) ".rst")
               out-root)))
        (make-directory (file-name-directory rst-file) t)
        (unless (= 0 (call-process "pandoc" nil "*cpmdc-pandoc*" t
                                   "-f" "org" "-t" "rst"
                                   "--lua-filter"
                                   (expand-file-name "org-links.lua"
                                                     default-directory)
                                   org-file "-o" rst-file))
          (with-current-buffer "*cpmdc-pandoc*"
            (princ (buffer-string)))
          (error "pandoc failed for %s" org-file))))))

(defun cpmdc-rst-file-doc-link (link desc info)
  "Export a relative Org file link as a Sphinx :doc: role.
INFO is accepted so this can advise `org-rst-link`. Return nil for every
other link so ox-rst handles it."
  (ignore info)
  (let ((type (org-element-property :type link))
        (path (org-element-property :path link)))
    (when (and (equal type "file")
               (stringp path)
               (string-match "\\`\\(.*\\)\\.org\\'" path))
      (let ((page (match-string 1 path))
            (text (if (and (stringp desc) (not (string-empty-p desc)))
                      desc
                    (file-name-nondirectory (match-string 1 path)))))
        (format ":doc:`%s <%s>`" text page)))))

(defun cpmdc-publish-with-ox-rst ()
  "Export Org sources with MELPA ox-rst."
  (require 'package)
  (add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
  (package-initialize)
  (unless (package-installed-p 'ox-rst)
    (package-refresh-contents)
    (package-install 'ox-rst))
  (require 'ox-rst)
  (require 'ox-publish)
  (advice-add 'org-rst-link :before-until #'cpmdc-rst-file-doc-link)
  (setq org-publish-project-alist
        '(("sphinx-rst"
           :base-directory "./orgmode/"
           :base-extension "org"
           :publishing-directory "./source/"
           :publishing-function org-rst-publish-to-rst
           :recursive t
           :headline-levels 4
           :with-author nil
           :with-toc nil
           :section-numbers nil)))
  (org-publish "sphinx-rst" t))

(if (string= (getenv "CPMDC_ORG_EXPORT") "pandoc")
    (cpmdc-publish-with-pandoc)
  (cpmdc-publish-with-ox-rst))

Application autonome Apple Silicon pour macOS 14 ou ultérieur : Python, noyau
Vinkulum 0.19.0, PySide6 6.11.2 et VTK 9.7.0 embarqués. Télécharger le DMG et
glisser Vinkulum Studio dans Applications. Aucune compilation locale nécessaire.

Version expérimentale, signature ad hoc, sans notarisation Apple. macOS peut
demander une autorisation d’ouverture dans Confidentialité et sécurité.

La construction vérifie les tests natifs et Studio sur ARM64, la signature,
l’architecture, le démarrage du bundle, le rendu OpenGL et un calcul du double
pendule par le worker embarqué. Le rapport et la capture sont joints. Ces
contrôles ne constituent pas une certification générale du noyau.

Le DMG est conservé comme asset de release, sans stockage volumineux dans les
artefacts GitHub Actions. La compilation utilise le runner public macos-14
standard. Sources exactes : tag de cette prerelease ; licences incluses au DMG.

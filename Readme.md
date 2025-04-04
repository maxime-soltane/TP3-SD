Vous commencez votre premier job de bioinformaticien.ne. La personne précédente vous a laissé un code quasi sans commentaire et sans test, et votre chef est biologiste et n'y comprend rien. Il sait seulement que ce code était très pratique pour rechercher des k-mer dans des jeux de données. 
Vous allez améliorer la situation pour que l'équipe puisse se servir du code de manière pérenne.

1. Décrivez à quoi sert `class SimpleBloomFilter`. Ecrivez les doctests pour cette classe.

La classe SimpleBloomFilter permet d'implémenter une structure appelée filtre de Bloom. Il s'agit d'une structure probabiliste permettant de vérifier si un élément, ici un kmer, est présent dans un jeu de données. Elle initialise un tableau de n bits, avec que des 0 et sa taille est définie par "size". Elle stocke les fonctions de hashage dont le nombre est défini par "num_hashes". Elle génère plusieurs positions de hashage via la methode "_hashes". Elle permet d'insérer un élément "item" dans le filtre de bloom par la méthode "add", pour cela, elle utilise le hash de l'item pour obtenir les positions et met à 1 toutes les positions correspondantes dans le "bit_array". Elle peut vérifier la probable présence d'un élément via la méthode "contains", elle renverra "True" si l'élément est probablement présent, avec des risques de faux positifs et "False" lorsque l'élément est certainement absent. Pour cela, elle recalcule les positions via la fonction "_hashes" et vérifie que toutes ses positions sont à 1. Grâce à la méthode "merge" on peut fusionner deux filtres de Bloom de même taille. Pour cela, la méthode compare les deux "bit_array" bit à bit et retourne un nouveau filtre, puis fait de l'union des deux filtres parents. Cette derniere méthode permet de construire une autre structure, un arbre de Bloom filters.

2. Décrivez à quoi sert `class StructureNode`.

La classe StructureNode modélise un noeud dans un arbre binaire construit à partir d'un ensemble de datasets, permettant d'indexer ces jeux de données via un filtre de Bloom (self.bloom). Le filtre de bloom permettra de vérifier la présence d'un kmer dans les datasets (self.datasets). La classe implémente également des objets "StructureNode" pour l'enfant de droite (self.right) et de gauche (self.left). Ainsi, chaque feuille de l'arbre contient un filtre de Bloom des kmers d'un dataset, et les noeuds parents sont produits par la fusion des filtres de bloom de leurs deux enfants, l'union de tous les kmers portés par l'ensemble des datasets correspond donc à la racine de l'arbre. Cette classe permet de faire une recherche dans l'arbre plus rapide et plus économique lorsque le kmer recherché n'est pas présent. Pour cela, l'algorithme part de la racine, vérifie la présence du kmer dans le filtre de bloom du noeud et descend l'arbre à gauche et à droite lorsqu'il est détecté. Si le noeud est une feuille, il retourne les datasets. Cela permet d'éviter de parcourir des branches entières si le filtre indique que le kmer est absent.

3. Décrivez à quoi sert `class Structure`. Ecrivez les doctests pour cette classe.

Cette classe permet d'organiser un ensemble de jeux de données (datasets) dans un arbre binaire. Chaque noeud est composé d'un filtre de bloom. Dans le constructeur on retrouve une liste de jeux de données (datasets), un dictionnaire pour chaque dataset et les paramètres de taille et de nombre de fonctions de hashage (bloom_size, num-hashes) pour la construction des filtres de bloom. La construction de l'arbre se fait grâce à la fonction "_build_tree()". Cette méthode crée une feuille par dataset, cette dernière sera composée d'un "StructureNode" et d'un filtre de bloom portant tous les kmers de ce dataset. Puis elle regroupe les noeuds en fusionnant leurs filtres via la méthode "merge", et ce jusqu'à arriver à la racine. De plus, la méthode "query" permet de vérifier dans quel dataset un kmer est probablement présent en parcourant l'arbre depuis la racine, puis noeud après noeud à gauche et à droite. Si le filtre contient le kmer, on continue, sinon on coupe dans cette branche. Lorsqu'une feuille est atteinte, on ajoute le nom du dataset à la réponse. L'emploi de filtres de Bloom permet de faire cela plus rapidement, sans parcourir entièrement les sous-arbres.

4. Cette structure mélange donc deux structures de données que nous avons vues. Quelles sont elles ?

Cette structure mélange les tables de hashage, puisque les filtres de Bloom sont des tables de hashage probabilistes, ainsi que la liste chaînée, dans sa structure. En effet, les listes chaînées sont faites d'éléments contenant une valeur et un lien vers l'élément suivant. Ici les objets StructureNode contiennent un filtre de Bloom et deux liens vers les noeuds suivants droit et gauche.

5. D'après vous, que peut-on dire sur la complexité de la requête de cette structure ? 

La complexité d'une requête dépend de la profodeur de l'arbre ainsi que du nombre de faux positifs, en principe, le parcours d'un arbre binaire équilibré se fait en O(log n). Cependant, lorsqu'on rencontre un grand nombre de faux positifs, l'algorithme peut alors parcourir beaucoup plus de branches inutilement, ce qui fait croître la compléxité. De plus, l'utilisation de filtres de Bloom permet d'éviter de très nombreuses vérifications dans le cas de très grand jeux de données, et limite ainsi la complexité.

6. Quelles sont les différences avec la table basée sur une MPHF que nous avons vu ?

Les différences entre ces deux structures sont :
- Les filtres de Bloom connaissent un risque de faux positifs contrairement aux MPHF qui sont exactes.
- L'espace utilisé par la MPHF est plus restreint que pour un filtre de Bloom.
- Un filtre de Bloom est malléable contrairement à une MPHF, on peut lui ajouter des données et les fusionner alors qu'une fois la MPHF construite, elle est statique.
- Les MPHF sont très compactes et rapides mais complexes à construire. Les filtres de Bloom nécessitet davantage d'espace, en raison des faux positifs, mais sont plus simples à construire et à maipuler.
- Les MPHF permettent de représenter des ensembles de kmers fixes et sans collision, alors que les filtres de Bloom sont plus flexibles et combinables et permettent de vérifier la présence d'un kmer dans des datasets probabilistiquement.

7. Bonus : Pouvez-vous retracer de quel papier de bioinformatique vient cette idée ?

# ICRH et très haute densité : CFS/SPARC et configuration type ITER

Nov 5, 2026 · @Robert Dosseh-Kpotogbey (Kpotogbey)

## 1. Thèse et périmètre

L'ICRH est la seule méthode de chauffage auxiliaire dont la pénétration ne se dégrade pas quand la densité monte, ce qui explique que CFS en fasse l'unique chauffage externe de SPARC. Le vrai verrou n'est pas la propagation dans le cœur, mais le couplage de l'antenne à travers le bord du plasma.

Cette note développe l'argument physique (section 2), le verrou du couplage (section 3), le choix de CFS (section 4), puis ce qu'exigerait une configuration ITER adaptée aux hautes densités (section 5).

**Lecture de votre question de précision.** Je l'interprète ainsi : quelles tolérances faut-il tenir (adaptation d'impédance, tension d'antenne, réponse aux ELMs, spectre lancé), et ITER, reconfiguré pour la haute densité, les tiendrait-il ? Si vous visiez autre chose, indiquez-le et j'ajusterai.

**Deux précisions sur mes réponses précédentes :**

- Pour SPARC, les publications de 2026 donnent 20 MW installés à 120 MHz pour la première campagne. Les études de 2018-2019 parlaient de 30 MW couplés. Mon chiffre d'environ 25 MW était donc approximatif.
- « Pas de densité de coupure » vaut pour la propagation dans le cœur. Il existe bien une couche évanescente au bord, là où la densité est trop faible pour propager l'onde rapide : c'est elle qui gouverne le couplage.

## 2. Physique de l'onde rapide : pourquoi la densité ne la bloque pas

L'antenne ICRF lance l'onde magnétosonore rapide, dont le champ électrique est surtout perpendiculaire au champ magnétique ([Lin et al., JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)). Pour un faible indice parallèle n∥, son nombre d'onde perpendiculaire suit approximativement la vitesse d'Alfvén (relation classique de plasma, à titre d'ordre de grandeur) :

```latex
k_\perp \simeq \frac{\omega}{v_A}, \qquad v_A = \frac{B}{\sqrt{\mu_0\, n\, m_i}}
```

Quand la densité n augmente, v\_A diminue et k⊥ croît comme √n : la longueur d'onde raccourcit, mais l'onde continue de se propager. Il n'existe pas de densité au-delà de laquelle elle serait réfléchie dans le cœur.

C'est le contraire des ondes électroniques, qui rencontrent une coupure dès que la fréquence plasma dépasse la fréquence de l'onde (section 6).

### Ce que disent les sources

Les concepteurs de SPARC l'écrivent explicitement : contrairement au NBI, le chauffage ICRF n'a pas de limite de densité pour la pénétration, et le point de fonctionnement à haute densité de SPARC favorise même le couplage ([Lin et al., JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)). Une présentation de l'école ITER 2026 résume l'ICRF ainsi : pas de limite de densité pour coupler la puissance, et une efficacité d'absorption qui augmente avec le bêta ionique ([Zhang, ITER International School 2026](https://www.iter.org/sites/default/files/media/2026-07/l-4_zhang.pdf)).

### La résonance fixe la fréquence, pas la densité

La fréquence cyclotronique d'une espèce est f = qB/(2πm). Elle ne dépend que du champ, pas de la densité. Pour SPARC (B₀ = 12,2 T), mon calcul donne environ 124 MHz pour l'hélium-3 fondamental sur l'axe, et la même valeur pour le tritium au 2ᵉ harmonique, car q/m de l'³He vaut presque le double de celui du T. Les deux absorptions se superposent donc, et une fréquence d'environ 120 MHz place la résonance près de l'axe.

Cette coïncidence est confirmée par les concepteurs : à ces champs, l'³He fondamental et le T au 2ᵉ harmonique résonnent vers 120 MHz ([Lin et al., JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)).

### Absorption en une passe

Un schéma de chauffage fonctionne bien s'il combine une forte absorption en simple passage (SPA) et une absorption localisée sur l'axe, d'après l'expérience d'Alcator C-Mod. Pour le plasma D-T de SPARC, les travaux de 2019 et 2020 annoncent une SPA supérieure à 50 % dans la plupart des conditions, et supérieure à 90 % en plateau ([APS 2019](https://archive.aps.org/dpp/2019/co5/8) ; [APS 2020](https://meetings-archive.aps.org/dpp/2020/jo08/11)). Ce sont des résultats de modélisation, pas des mesures.

Les simulations 2026 pour la première campagne montrent que la SPA reste peu sensible à la concentration d'³He entre 3 % et 5 %, et qu'elle baisse quand l'³He diminue ([Migliore et al., EPJ Web Conf. 2026](https://www.epj-conferences.org/10.1051/epjconf/202634602011)). Autrement dit, la haute densité n'est pas le facteur limitant ; le choix de la concentration de minoritaire l'est davantage.

## 3. Le vrai verrou : le couplage antenne-plasma

La limite de l'ICRH est la couche évanescente entre l'antenne et la coupure de l'onde rapide au bord du plasma. La résistance de charge de l'antenne dépend de la distance entre l'antenne et cette coupure, donc du profil de densité dans la zone de bord ([Mayoral et al., IAEA 2010](https://www-pub.iaea.org/MTCD/Meetings/PDFplus/2010/cn180/cn180_papers/itr_p1-11.pdf) ; [Messiaen et al.](https://www.ipfn.tecnico.ulisboa.pt/soft2010_abstractbook/pdf/P2-112.pdf)).

La puissance couplée est de la forme suivante, avec R\_c la résistance de charge et I le courant dans la bande :

```latex
P \simeq \tfrac{1}{2}\, R_c\, I^2
```

À puissance égale, une charge faible impose un courant et une tension plus élevés. On se rapproche alors des limites de claquage de l'antenne. Cette relation est générale ; je ne l'ai pas tirée d'un des documents cités.

### Les ELMs fixent la limite réelle

En mode H, ce n'est pas la charge moyenne qui limite la puissance, mais la charge minimale entre deux ELMs ([EPS 2010, O4.124](https://pure.mpg.de/pubman/item/item_2140803_3/component/file_3325673/O4.124.pdf)). Sur JET, une injection locale de gaz a fait passer cette charge minimale d'environ 0,17 Ω à environ 1 Ω. Sur DIII-D, le gain a atteint un facteur 6 ([Mayoral et al.](https://www-pub.iaea.org/MTCD/Meetings/PDFplus/2010/cn180/cn180_papers/itr_p1-11.pdf)).

Le même article EPS conclut que l'incertitude sur la charge attendue dans ITER vient de la grande variété possible de profils de densité dans la SOL lointaine. Les modèles 3D comme TOPICA reproduisent bien les mesures quand le profil local est connu.

### Ce qui change à haute densité

Une densité de bord plus élevée rapproche la coupure de l'antenne et réduit l'évanescence. Les concepteurs de SPARC l'affirment : la charge d'antenne de SPARC, fixée surtout par le profil de densité au bord, est proche de celle d'Alcator C-Mod et bien supérieure à celle d'ITER ([Lin, APS 2018](https://archive.aps.org/dpp/2018/go5/4)). Leur résumé parle d'une charge modérée pour un bon couplage, associée à une forte absorption en simple passage.

La haute densité aide donc le couplage au lieu de le pénaliser. Elle déplace toutefois la difficulté vers l'ingénierie du bord.

### Trois contraintes qui restent

- **Variations de charge.** SPARC prévoit des tuners à ferrite pour l'adaptation en temps réel ([Lin, APS 2018](https://archive.aps.org/dpp/2018/go5/4)), et un contrôleur qui ajuste l'adaptation en cours d'impulsion par modulation de fréquence de ±1 MHz ([Usoltseva et al., EPJ Web Conf. 2026](https://doi.org/10.1051/epjconf/202634602001)).
- **Impuretés.** Dans une machine toute en métal, la pulvérisation du tungstène peut limiter la puissance ICRF utilisable. Les mesures de réduction passent par l'alignement des antennes sur le champ, la minimisation des courants image et la modélisation des gaines RF ([Lin, APS 2018](https://archive.aps.org/dpp/2018/go5/4) ; [UKAEA, 2026](https://indico.ukaea.uk/event/449/contributions/1743/)).
- **Spectre lancé.** Le spectre de k∥ de l'antenne influence à la fois le couplage et l'absorption ; il se règle par conception et doit être optimisé pour chaque machine.

## 4. CFS/SPARC : un tokamak construit autour de l'ICRF

SPARC est un tokamak compact à haut champ et haute densité, dont le seul chauffage externe est l'ICRF. Les paramètres de conception publiés sont B₀ = 12,2 T, I\_p = 8,7 MA, R = 1,85 m, a = 0,57 m, n\_e0 ≈ 4×10²⁰ m⁻³ et T\_e0 ≈ 20 keV ([Lin et al., JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)).

### Pourquoi pas le NBI, ni l'ECRH

Pour pénétrer le plasma dense de SPARC, les ions neutres devraient atteindre des énergies du niveau du MeV, peut-être au-delà du faisceau prévu pour ITER ([Lin et al., JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)). L'article ajoute que les autres méthodes présentent toutes des difficultés importantes de physique, de coût ou de technologie. Il retient en plus la disponibilité de sources RF à haut champ et l'absence de limite de densité pour la pénétration.

La taille de la machine compte aussi : les ports sont limités, et l'ICRF se loge dans des antennes montées sur la paroi interne du vaisseau.

### Schémas de chauffage et fréquence

L'ICRF de SPARC fonctionne à 120 ± 1 MHz. Deux scénarios sont étudiés pour la première campagne ([UKAEA, 2026](https://indico.ukaea.uk/event/449/contributions/1743/)) :

- **Environ 12 T :** minoritaire ³He, avec absorption harmonique sur le tritium.
- **Environ 8,5 à 8,9 T :** minoritaire H, avec absorption harmonique sur D et T.

Les simulations 2026 valident le schéma ³He + T au 2ᵉ harmonique pour le scénario L-mode Q > 1, avec un chauffage dominé par les ions et une absorption efficace en simple passage ([Usoltseva et al., EPJ Web Conf. 2026](https://doi.org/10.1051/epjconf/202634602001)).

### Matériel : ce qui est annoncé

| Élément | Valeur annoncée | Source |
| --- | --- | --- |
| Puissance installée, 1ʳᵉ campagne | 20 MW à 120 MHz | [Usoltseva et al. 2026](https://doi.org/10.1051/epjconf/202634602001) |
| Antennes (présentation d'ensemble) | 10 antennes à 4 bandes, en paires poloïdales | [Usoltseva et al. 2026](https://doi.org/10.1051/epjconf/202634602001) |
| Antennes (simulation 2026) | 14 antennes à 7 positions toroïdales, plus de 20 MW | [EPJ Web Conf. 02011](https://www.epj-conferences.org/10.1051/epjconf/202634602011) |
| Plan de 2020 | 12 antennes à 4 bandes dans 6 ports, avec antennes de secours à 3 bandes | [Lin et al., JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader) |
| Réglage en impulsion | Modulation de fréquence de ±1 MHz | [Usoltseva et al. 2026](https://doi.org/10.1051/epjconf/202634602001) |

Les nombres d'antennes diffèrent d'un document à l'autre. Il s'agit probablement d'évolutions de conception, mais je n'ai pas de source qui l'établisse.

### Où en est SPARC

Selon les documents consultés, SPARC est en assemblage et en préparation de sa première campagne. Toutes les valeurs d'absorption, de charge et de puissance ci-dessus sont des valeurs de conception ou de modélisation, validées par l'expérience d'Alcator C-Mod et de JET, pas par des plasmas SPARC ([UKAEA, 2026](https://indico.ukaea.uk/event/449/contributions/1743/)). Je n'ai pas vérifié si des plasmas ont été produits depuis.

## 5. ITER : la précision requise, et une configuration adaptée aux hautes densités

ITER et SPARC ont le même principe de chauffage, mais leurs difficultés de couplage sont opposées. SPARC a une charge d'antenne élevée et un risque d'impuretés ; ITER a une charge faible et une marge de tension étroite.

### Le système ICRH d'ITER

Le système IC H&CD d'ITER doit coupler 20 MW en régime quasi continu (impulsions jusqu'à 3600 s) entre 40 et 55 MHz, avec un couplage robuste pendant les ELMs. Il comprend deux antennes en ports équatoriaux, chacune dimensionnée pour 20 MW, avec une limite de tension de 45 kV ([AIP Conf. Proc. 1187, p. 265](https://scientific-publications.ukaea.uk/wp-content/uploads/Published/AIPCONFPROC1187p265.pdf)). Chaque antenne est un réseau de 24 bandes (4 toroïdales × 6 poloïdales) ([Goulding & Swain, APS 2006](https://archive.aps.org/dpp/2006/vp1/50)). Le système est extensible jusqu'à 40 MW ([Fusion Eng. Des. 2013](https://bvbvf1.bib-bvb.de/vufind/PrimoRecord/cdi_hal_primary_oai_HAL_cea_04748556v1)).

À 5,3 T, le schéma tritium 2ᵉ harmonique + minoritaire He se situe vers 53-55 MHz selon le tableau de scénarios de la même publication. Mon calcul donne environ 54 MHz pour le T au 2ᵉ harmonique, cohérent avec ce tableau.

### SPARC et ITER côte à côte

| Grandeur | SPARC | ITER |
| --- | --- | --- |
| Champ toroïdal | 12,2 T (8,5 à 8,9 T en champ réduit) | 5,3 T nominal |
| Densité centrale | ≈ 4×10²⁰ m⁻³ | ≈ 1×10²⁰ m⁻³ (ordre de grandeur, non sourcé ici) |
| Fréquence ICRF | 120 ± 1 MHz | 40 à 55 MHz |
| Puissance ICRF | 20 MW installés (1ʳᵉ campagne) | 20 MW, extensible à 40 MW |
| Rôle | Seul chauffage externe | Un des trois systèmes, avec NBI et ECRH |
| Antennes | 10 ou 14 selon la source, 4 bandes | 2 antennes, 24 bandes chacune |
| Charge d'antenne | Proche de C-Mod, bien plus haute qu'ITER | Faible : séparatrice à environ 14 cm de l'antenne |
| Tenue aux ELMs | Tuners à ferrite, ±1 MHz de modulation | Réseau passif résilient à la charge, 45 kV max |

Sources : [Lin 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader), [Usoltseva 2026](https://doi.org/10.1051/epjconf/202634602001), [Lin APS 2018](https://archive.aps.org/dpp/2018/go5/4), [ITER](https://www.iter.org/node/15918), [AIP 1187](https://scientific-publications.ukaea.uk/wp-content/uploads/Published/AIPCONFPROC1187p265.pdf), [JET/ITER, OSTI](https://www.osti.gov/biblio/21035923).

L'absorption en simple passage de SPARC est proche de celle d'ITER, car la taille de la machine la renforce ; elle est plus élevée que celle de C-Mod ([Lin, APS 2018](https://archive.aps.org/dpp/2018/go5/4)). Sur ce critère, un ITER dense et un SPARC se ressemblent.

### Les tolérances à tenir côté ITER

- **Source.** Une unité doit fournir 2,5 MW sur un ROS (VSWR) de 2, pour toute phase du coefficient de réflexion, et 3 MW sur un ROS de 1,5 dans la bande de fonctionnement ([Kazarian et al., SOFT 2010](https://www.ipfn.tecnico.ulisboa.pt/soft2010_abstractbook/pdf/P1-093.pdf)).
- **Charge.** Dans ITER, l'antenne doit délivrer 20 MW en mode H avec ELMs, pour une distance moyenne antenne-séparatrice de 14 cm. Deux problèmes se posent : la variation très rapide de la charge pendant les ELMs et sa chute quand le plasma s'éloigne ([OSTI](https://www.osti.gov/biblio/21035923)).
- **Tension.** La grande distance implique une faible résistance de charge, donc une tenue en tension proche du haut de ce qui a été démontré ([ITER, couplage](https://core.ac.uk/works/4833398)).

### Une configuration ITER adaptée à la haute densité : mon analyse

ITER n'a pas, à ma connaissance, de programme officiel de ce type. Ce qui suit est mon raisonnement à partir des sources ci-dessus, pas une position de l'Organisation ITER.

1. **Fréquence et schéma inchangés.** À 5,3 T, la bande de 40 à 55 MHz couvre déjà le schéma T 2ᵉ harmonique. La densité ne déplace pas la résonance ; elle renforce plutôt l'absorption par le bêta ionique.
2. **Charge en hausse.** Une densité de séparatrice plus élevée, obtenue par injection de gaz près de l'antenne, a multiplié la charge entre ELMs par environ 6 sur JET et jusqu'à 6 sur DIII-D ([Mayoral et al.](https://www-pub.iaea.org/MTCD/Meetings/PDFplus/2010/cn180/cn180_papers/itr_p1-11.pdf)). Cela corrige la faiblesse principale d'ITER. Densité du cœur et densité de bord ne sont toutefois pas identiques : un plasma dense au centre n'a pas forcément un bord dense.
3. **Marge de tension.** Avec P ∝ R\_c I², une charge plus forte réduit le courant et la tension pour 20 MW. Cela donnerait de la marge sous la limite de 45 kV et pourrait aider une montée vers 40 MW. Le refroidissement et les alimentations restent des limites que je n'ai pas évaluées.
4. **Adaptation active.** ITER utilise un réseau passif résilient à la charge. SPARC ajoute un réglage actif en fréquence de ±1 MHz. Cette option serait techniquement concevable dans une bande aussi large que 40-55 MHz, mais je n'ai pas de source indiquant qu'elle soit envisagée.
5. **Coût en confinement.** Sur JET, le puff de gaz a augmenté la fréquence des ELMs (jusqu'à environ 100 Hz) et dégradé le confinement global ([EPS 2010, O4.124](https://pure.mpg.de/pubman/item/item_2140803_3/component/file_3325673/O4.124.pdf)). La réponse n'est donc pas gratuite : on échange du couplage contre de la performance.

### Réponse à votre question

La précision requise est de même nature chez ITER et chez SPARC (tenue au ROS, aux ELMs et à la tension), mais pas de même valeur. Chez ITER elle est fixée par une charge faible et une tension limitée. Chez SPARC elle l'est par une charge forte, la modulation de fréquence et la gestion du tungstène. Un ITER configuré pour la haute densité se rapprocherait du régime de SPARC sur la charge, sans le rejoindre sur la taille de la machine.

## 6. Comparaison à haute densité : ICRH, ECRH, LHCD, NBI

À haute densité, l'ICRH est la seule méthode dont aucune limite physique de pénétration n'est atteinte ; les trois autres rencontrent chacune un seuil ou une perte d'efficacité.

| Méthode | Ce qui limite à haute densité | Position de SPARC et d'ITER |
| --- | --- | --- |
| ICRH/ICRF | Le couplage de bord, pas la propagation : aucune limite de densité pour la pénétration ([Lin 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)) | SPARC : seul chauffage externe. ITER : 20 MW, 40-55 MHz |
| ECRH/ECH | Coupure de l'onde ordinaire : n\_c ≈ 1,24×10⁻² f² (m⁻³, f en Hz), soit environ 2,4×10²⁰ m⁻³ à 140 GHz. Atteindre 4×10²⁰ m⁻³ exigerait environ 180 GHz (mon calcul) | SPARC : l'article de conception indique des difficultés de physique, de coût ou de technologie. ITER : 170 GHz, 20 MW (valeurs de ma réponse précédente, non revérifiées) |
| LHCD | Absorption dans le bord à haute densité, donc perte de génération de courant (connaissance générale, non sourcée ici) | Absent d'ITER |
| NBI | Les neutres s'ionisent en périphérie : il faut des énergies de l'ordre du MeV. Pour SPARC, peut-être au-delà du faisceau d'ITER ([Lin 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)) | ITER : 1 MeV D⁻, 33 MW (d'après votre texte) |

L'ordre de préférence à très haute densité est donc : ICRH, puis ECRH avec ses limites de coupure, puis NBI à très haute énergie, puis LHCD. C'est celui que j'avais donné précédemment ; les sources de SPARC le confirment pour les deux extrêmes (ICRH et NBI), et le calcul de coupure le confirme pour l'ECRH.

## 7. Limites et points à vérifier

- **Valeurs de conception, pas mesures.** Pour SPARC, absorption, charge d'antenne et puissance viennent de la modélisation et de l'expérience de C-Mod et JET. Aucune donnée de plasma SPARC n'apparaît dans les sources consultées.
- **Chiffres qui varient.** La puissance est passée de 30 MW couplés (2018-2019) à 20 MW installés (2026). Le nombre d'antennes est de 10, 12 ou 14 selon le document.
- **Sources ITER anciennes.** Plusieurs références ITER datent de 2006 à 2017. Je n'ai pas vérifié l'état actuel de la conception ICRH ni le calendrier de ses phases.
- **Calculs personnels.** Les fréquences de résonance (124 MHz pour SPARC, environ 54 MHz pour ITER) et les densités de coupure ECRH sont mes calculs, signalés comme tels dans le texte. Les rangs de l'ECRH et du LHCD reposent en partie sur ma connaissance générale.
- **Absorption non garantie.** L'absence de limite de densité n'assure pas l'absorption : celle-ci dépend de la concentration de minoritaire, de la température et du bêta ionique, et une part de la puissance peut chauffer les électrons.
- **Densité de bord.** Une haute densité de cœur ne fixe pas la densité de séparatrice. Le gain de charge par injection de gaz a un coût en confinement, mesuré sur JET.
- **Consultation des sources.** J'ai travaillé à partir d'extraits de pages et de résumés. Je n'ai pas ouvert chaque article en entier, donc les valeurs de détail sont à confirmer avant publication.

## 8. Conclusion

Avec CFS, l'ICRH est le chauffage qui s'accommode le mieux de la très haute densité, parce que l'onde rapide se propage dans le cœur à toute densité et que l'absorption croît avec le bêta ionique. Sa limite est le couplage de bord, que la haute densité allège plutôt qu'elle ne l'aggrave.

Pour ITER, une configuration adaptée à la haute densité garderait la bande de 40 à 55 MHz et le schéma T 2ᵉ harmonique. Elle chercherait à relever la charge d'antenne par le profil de densité de bord, au prix d'une dégradation du confinement. Les tolérances à tenir (ROS, ELMs, tension) restent de même nature que chez SPARC, avec des valeurs fixées par une autre physique de bord.

La validation décisive viendra des premières campagnes de SPARC, qui mesureront enfin charge, absorption et impuretés en régime de très haute densité.

## Sources

- [Physics basis for the ICRF system of the SPARC tokamak, JPP 2020](https://core-prod.cambridgecore.org/core/product/22016DD64F3C5CAD47563A1E4AE59934/core-reader)
- [Overview of the ICRF heating system in SPARC, EPJ Web Conf. 2026](https://doi.org/10.1051/epjconf/202634602001)
- [Full-wave ICRF modelling for SPARC, EPJ Web Conf. 2026](https://www.epj-conferences.org/10.1051/epjconf/202634602011)
- [ICRF power coupling and absorption in SPARC's first campaign, UKAEA 2026](https://indico.ukaea.uk/event/449/contributions/1743/)
- [ICRF Heating for SPARC, APS DPP 2018](https://archive.aps.org/dpp/2018/go5/4) ; [2019](https://archive.aps.org/dpp/2019/co5/8) ; [2020](https://meetings-archive.aps.org/dpp/2020/jo08/11)
- [Zhang, ITER International School 2026](https://www.iter.org/sites/default/files/media/2026-07/l-4_zhang.pdf)
- [ITER ion cyclotron heating](https://www.iter.org/node/15918)
- [ITER IC H&CD system requirements, AIP Conf. Proc. 1187](https://scientific-publications.ukaea.uk/wp-content/uploads/Published/AIPCONFPROC1187p265.pdf)
- [Advances in ITER Ion Cyclotron System Design, APS 2006](https://archive.aps.org/dpp/2006/vp1/50)
- [Status of the ITER IC H&CD system, Fusion Eng. Des. 2013](https://bvbvf1.bib-bvb.de/vufind/PrimoRecord/cdi_hal_primary_oai_HAL_cea_04748556v1)
- [RF sources for ITER IC H&CD, SOFT 2010](https://www.ipfn.tecnico.ulisboa.pt/soft2010_abstractbook/pdf/P1-093.pdf)
- [On maximizing the ICRF antenna loading for ITER plasmas, IAEA 2010](https://www-pub.iaea.org/MTCD/Meetings/PDFplus/2010/cn180/cn180_papers/itr_p1-11.pdf)
- [Gas puffing and ICRF coupling, EPS 2010](https://pure.mpg.de/pubman/item/item_2140803_3/component/file_3325673/O4.124.pdf)
- [Simulation of ICRF antenna loading by a dielectric dummy load, SOFT 2010](https://www.ipfn.tecnico.ulisboa.pt/soft2010_abstractbook/pdf/P2-112.pdf)
- [Coupling of the JET ICRF antennas with ITER-relevant plasma-strap distance, OSTI](https://www.osti.gov/biblio/21035923)

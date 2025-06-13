import arcadia_pycolor as apc

PLOTLY_TITLE_FONT = apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold"

GRAY_GRADIENT = (
    apc.Gradient(name="grays", colors=[apc.slate, apc.charcoal, apc.steel, apc.cloud, apc.dove])
    .interpolate_lightness()
    .reverse()
)

CINNABAR_GRADIENT = apc.Gradient(
    name="cinnabars", colors=[apc.melon, apc.dragon, apc.cinnabar, apc.redwood]
).interpolate_lightness()


DEFAULT_FUNNEL_COLOR_LIST = [
    apc.aster,
    apc.lapis,
    apc.aegean,
    apc.vital,
    apc.seaweed,
    apc.teal,
    apc.glass,
    apc.matcha,
    apc.lime,
]

aegeanlike = apc.HexCode("aegeanlike", "#2C65A1")
iris = apc.HexCode("iris", "#DCBFFC")

DEFAULT_SPECIES_COLORS = {
    "Abeoforma-whisleri": apc.mars,
    "Aedes-aegypti": apc.dusk,
    "Agaricus-bisporus": apc.taupe,
    "Amorphochlora-amoebiformis": apc.glass,
    "Anolis-carolinensis": apc.tanzanite,
    "Aspergillus-nidulans": apc.denim,
    "Bathycoccus-prasinos": apc.vital,
    "Bodo-saltans": apc.sky,
    "Caenorhabditis-elegans": apc.sun,
    "Callithrix-jacchus": apc.tanzanite,
    "Callorhinchus-milii": apc.tanzanite,
    "Candida-albicans": apc.denim,
    "Carlito-syrichta": apc.tanzanite,
    "Chlamydomonas-reinhardtii": apc.dragon,
    "Chlorella-vulgaris": apc.chateau,
    "Ciona-intestinalis": apc.iris,
    "Clytia-hemisphaerica": apc.cinnabar,
    "Danio-rerio": apc.tanzanite,
    "Dictyostelium-discoideum": apc.rose,
    "Diplonema-papillatum": apc.marine,
    "Drosophila-melanogaster": apc.dusk,
    "Entamoeba-histolytica": apc.aster,
    "Euglena-gracilis": apc.tangerine,
    "Exaiptasia-diaphana": apc.cinnabar,
    "Gallus-gallus": apc.tanzanite,
    "Giardia-intestinalis": apc.dress,
    "Hofstenia-miamia": aegeanlike,
    "Homo-sapiens": apc.tanzanite,
    "Hydra-vulgaris": apc.cinnabar,
    "Isochrysis-galbana": apc.oat,
    "Macaca-mulatta": apc.tanzanite,
    "Microcebus-murinus": apc.tanzanite,
    "Micromonas-commoda": apc.vital,
    "Mnemiopsis-leidyi": apc.melon,
    "Monosiga-brevicollis": iris,
    "Mus-musculus": apc.tanzanite,
    "Naegleria-gruberi": apc.sage,
    "Nannochloropsis-sp": apc.asparagus,
    "Nematostella-vectensis": apc.cinnabar,
    "Neurospora-crassa": apc.denim,
    "Ostreococcus-tauri": apc.vital,
    "Pan-troglodytes": apc.tanzanite,
    "Paramecium-tetraurelia": apc.amber,
    "Penicillium-chrysogenum": apc.denim,
    "Perkinsus-marinus": apc.lime,
    "Petromyzon-marinus": apc.tanzanite,
    "Phaeodactylum-tricornutum": apc.asparagus,
    "Plasmodium-falciparum": apc.aegean,
    "Porphyra-yezoensis": apc.teal,
    "Pristionchus-pacificus": apc.sun,
    "Saccharomyces-cerevisiae": apc.denim,
    "Salpingoeca-rosetta": iris,
    "Schistosoma-mansoni": apc.umber,
    "Schizosaccharomyces-pombe": apc.denim,
    "Schmidtea-mediterranea": apc.umber,
    "Sphaeroforma-arctica": apc.mars,
    "Symbiodinium-sp": apc.canary,
    "Taeniopygia-guttata": apc.tanzanite,
    "Tetrahymena-thermophila": apc.amber,
    "Tetraselmis-striata": apc.seaweed,
    "Ustilago-maydis": apc.taupe,
    "Volvox-carteri": apc.dragon,
    "Xenopus-tropicalis": apc.tanzanite,
    "Yarrowia-lipolytica": apc.denim,
}

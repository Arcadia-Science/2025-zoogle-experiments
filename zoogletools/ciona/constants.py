from enum import StrEnum

import arcadia_pycolor as apc


# The correct order of developmental stages is as follows:
class CionaStage(StrEnum):
    INIG = "iniG"
    MIDG = "midG"
    EARN = "earN"
    LATN = "latN"
    INITI = "iniTI"
    EARTI = "earTI"
    MIDTII = "midTII"
    LATTI = "latTI"
    LATTII = "latTII"
    LARVA = "larva"

    @classmethod
    def ordered_stages(cls):
        return [
            cls.INIG,
            cls.MIDG,
            cls.EARN,
            cls.LATN,
            cls.INITI,
            cls.EARTI,
            cls.MIDTII,
            cls.LATTI,
            cls.LATTII,
            cls.LARVA,
        ]


# Some of the stages are not named correctly in the Piekarz repo.
# The key is the Cao stage name, and the value is the (sometimes incorrect) Piekarz stage name.
# For example, if you want "latN" from Cao, you have to use "iniTI" file from Piekarz.
CIONA_STAGE_CAO_TO_PIEKARZ_MAP = {
    CionaStage.INIG: CionaStage.INIG,
    CionaStage.MIDG: CionaStage.MIDG,
    CionaStage.EARN: CionaStage.EARN,
    CionaStage.LATN: CionaStage.INITI,
    CionaStage.INITI: CionaStage.EARTI,
    CionaStage.EARTI: CionaStage.LATN,
    CionaStage.MIDTII: CionaStage.MIDTII,
    CionaStage.LATTI: CionaStage.LATTI,
    CionaStage.LATTII: CionaStage.LATTII,
    CionaStage.LARVA: CionaStage.LARVA,
}
CIONA_STAGE_PIEKARZ_TO_CAO_MAP = {v: k for k, v in CIONA_STAGE_CAO_TO_PIEKARZ_MAP.items()}

# This is a mapping between the notation used in the Cao et al. paper for developmental stages
# and the notation used in the Piekarz et al. paper for developmental stages.
CAO_STAGE_NAME_TO_PIEKARZ_STAGE_NAME_MAP = {
    "C110": CionaStage.INIG,
    "midG": CionaStage.MIDG,
    "earlyN": CionaStage.EARN,
    "lateN": CionaStage.LATN,
    "ITB": CionaStage.INITI,
    "ETB": CionaStage.EARTI,
    "MTB": CionaStage.MIDTII,
    "LTB1": CionaStage.LATTI,
    "LTB2": CionaStage.LATTII,
    "lv": CionaStage.LARVA,
}

# Colors used for developmental stages in the disambiguation notebook.
STAGE_COLORS = {
    CionaStage.INIG: apc.rose,
    CionaStage.MIDG: apc.dragon,
    CionaStage.EARN: apc.amber,
    CionaStage.LATN: apc.tangerine,
    CionaStage.INITI: apc.canary,
    CionaStage.EARTI: apc.lime,
    CionaStage.MIDTII: apc.seaweed,
    CionaStage.LATTI: apc.vital,
    CionaStage.LATTII: apc.aegean,
    CionaStage.LARVA: apc.aster,
}


def smooth_gradient_from_palette(palette: apc.Palette):
    return apc.Gradient(palette.name, palette.colors).interpolate_lightness()


PALETTE_DICT = {
    "epidermis": smooth_gradient_from_palette(apc.palettes.blue_shades),
    "nervous_system": smooth_gradient_from_palette(apc.palettes.purple_shades),
    "notochord": smooth_gradient_from_palette(apc.palettes.warm_gray_shades),
    "mesenchyme": smooth_gradient_from_palette(apc.palettes.red_shades),
    "muscle-heart": smooth_gradient_from_palette(apc.palettes.pink_shades),
    "endoderm": smooth_gradient_from_palette(apc.palettes.yellow_shades),
    "unannotated": smooth_gradient_from_palette(apc.palettes.cool_gray_shades),
    "germ": smooth_gradient_from_palette(apc.palettes.green_shades),
}

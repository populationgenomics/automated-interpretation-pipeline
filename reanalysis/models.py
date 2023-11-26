"""
A home for all data models used in AIP
"""


from pydantic import BaseModel, Field
from reanalysis.utils import get_granular_date


NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM


class Coordinates(BaseModel):
    """
    A representation of genomic coordinates
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    def string_format(self) -> str:
        """
        forms a string representation: chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def __lt__(self, other) -> bool:
        """
        enables positional sorting
        """
        # this will return False for same chrom and position
        if self.chrom == other.chrom:
            return self.pos < other.pos
        # otherwise take the relative index from sorted chromosomes list
        if self.chrom in CHROM_ORDER and other.chrom in CHROM_ORDER:
            return CHROM_ORDER.index(self.chrom) < CHROM_ORDER.index(other.chrom)
        # if self is on a canonical chromosome, sort before HLA/Decoy etc.
        if self.chrom in CHROM_ORDER:
            return True
        return False

    def __eq__(self, other) -> bool:
        """
        equivalence check
        Args:
            other (Coordinates):

        Returns:
            true if self == other

        """
        return (
            self.chrom == other.chrom
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
        )


class Variant(BaseModel):
    """
    the abstracted representation of a variant from any source

    Discriminators could be interesting to separate SV/Small/STR
    https://docs.pydantic.dev/latest/concepts/fields/#discriminator
    """

    info: dict[str, str | int | float]
    coordinates: Coordinates = Field(repr=True)
    categories: list[str] = Field(default_factory=list)
    het_samples: set[str] = Field(default_factory=set, exclude=True)
    hom_samples: set[str] = Field(default_factory=set, exclude=True)
    boolean_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_categories: list[str] = Field(default_factory=list, exclude=True)
    sample_support: list[str] = Field(default_factory=list, exclude=True)
    phased: dict = Field(default_factory=dict)

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coords < other.coords

    def __eq__(self, other):
        return self.coords == other.coords

    @property
    def has_boolean_categories(self) -> bool:
        """
        check that the variant has at least one assigned class
        """
        return any(self.info[value] for value in self.boolean_categories)

    @property
    def has_sample_categories(self) -> bool:
        """
        check that the variant has any list-category entries
        """
        return any(self.info[value] for value in self.sample_categories)

    @property
    def has_support(self) -> bool:
        """
        check for a True flag in any CategorySupport* attribute
        Returns:
            True if variant is support
        """
        return any(self.info[value] for value in self.sample_support)

    @property
    def category_non_support(self) -> bool:
        """
        check the variant has at least one non-support category assigned
        Returns:
            True if var has a non-support category assigned
        """
        return self.has_sample_categories or self.has_boolean_categories

    @property
    def is_classified(self) -> bool:
        """
        check for at least one assigned class, inc. support
        Returns:
            True if classified
        """
        return self.category_non_support or self.has_support

    @property
    def support_only(self) -> bool:
        """
        check that the variant is exclusively cat. support
        Returns:
            True if support only
        """
        return self.has_support and not self.category_non_support

    def sample_support_only(self, sample_id: str) -> bool:
        """
        check that the variant is exclusively cat. support
        check that this sample is missing from sample flags

        Returns:
            True if support only
        """
        return self.has_support and not self.sample_categorised_check(sample_id)

    def category_values(self, sample: str) -> list[str]:
        """
        get all variant categories
        steps category flags down to booleans - true for this sample

        Args:
            sample (str): sample id

        Returns:
            list of all categories applied to this variant
        """

        # step down all category flags to boolean flags
        categories = [
            category.replace('categorysample', '')
            for category in self.sample_categories
            if sample in self.info[category]  # type: ignore
        ]
        categories.extend(
            [
                bool_cat.replace('categoryboolean', '')
                for bool_cat in self.boolean_categories
                if self.info[bool_cat]
            ]
        )

        if self.has_support:
            categories.append('support')

        return categories

    def sample_categorised_check(self, sample_id: str) -> bool:
        """
        check if any *sample categories applied for this sample

        Args:
            sample_id (str): the specific sample ID to check

        Returns:
            bool: True if this sample features in any
                  named-sample category, includes 'all'
        """

        return any(
            sam in self.info[sam_cat]  # type: ignore
            for sam_cat in self.sample_categories
            for sam in [sample_id, 'all']
        )

    def sample_category_check(self, sample_id: str, allow_support: bool = True) -> bool:
        """
        take a specific sample and check for assigned categories
        optionally, include checks for support category

        Args:
            sample_id (str):
            allow_support: (bool) also check for support

        Returns:
            True if the variant is categorised for this sample
        """
        big_cat = self.has_boolean_categories or self.sample_categorised_check(
            sample_id
        )
        if allow_support:
            return big_cat or self.has_support
        return big_cat


class SmallVariant(Variant):
    depths: dict[str, int] = Field(default_factory=dict, exclude=True)
    ab_ratios: dict[str, float] = Field(default_factory=dict, exclude=True)
    transcript_consequences: list[dict[str, str]] = Field(default_factory=list)

    def get_sample_flags(self, sample: str) -> list[str]:
        """
        gets all report flags for this sample - currently only one flag
        """
        return self.check_ab_ratio(sample) + self.check_read_depth(sample)

    def check_read_depth(self, sample: str, threshold: int = 10) -> list[str]:
        """
        flag low read depth for this sample

        Args:
            sample (str): sample ID matching VCF
            threshold (int): cut-off for flagging as a failure

        Returns:
            return a flag if this sample has low read depth
        """
        if self.depths[sample] < threshold:
            return ['Low Read Depth']
        return []

    def check_ab_ratio(self, sample: str) -> list[str]:
        """
        AB ratio test for this sample's variant call

        Args:
            sample (str): sample ID

        Returns:
            list[str]: empty, or indicating an AB ratio failure
        """
        het = sample in self.het_samples
        hom = sample in self.hom_samples
        variant_ab = self.ab_ratios.get(sample, 0.0)
        if (
            (variant_ab <= 0.15)
            or (het and not 0.25 <= variant_ab <= 0.75)
            or (hom and variant_ab <= 0.85)
        ):
            return ['AB Ratio']
        return []


class SV_Variant(Variant):
    def check_ab_ratio(self, *args, **kwargs) -> list[str]:
        """
        dummy method for AB ratio checking - not implemented for SVs
        """
        return []

    def get_sample_flags(self, *args, **kwargs) -> list[str]:
        """
        dummy method for flag checking - not implemented for SVs (yet)
        """
        return []

    def check_read_depth(self, *args, **kwargs) -> list[str]:
        """
        dummy method for read depth checking - not implemented for SVs
        """
        return []


class ReportVariant(BaseModel):
    """
    A representation of a variant in a report
    """

    var_data: SmallVariant | SV_Variant
    sample: str
    categories: list[str] = Field(default_factory=list)
    family: str = Field(default_factory=str)
    gene: str = Field(default_factory=str)
    reasons: set[str] = Field(default_factory=set)
    genotypes: dict[str, str] = Field(default_factory=dict)
    support_vars: set[str] = Field(default_factory=set)
    flags: list[str] = Field(default_factory=list)
    panels: dict[str, str | list[int]] = Field(default_factory=dict)
    phenotypes: list[str] = Field(default_factory=list)
    labels: list[str] = Field(default_factory=list)
    first_seen: str = Field(default=get_granular_date())
    independent: bool = False

    @property
    def is_independent(self):
        """
        check if this variant acts independently
        """
        return len(self.support_vars) == 0

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        # self_supvar = set(self.support_vars)
        # other_supvar = set(other.support_vars)
        return (
            self.sample == other.sample
            and self.var_data.coords == other.var_data.coords
        )

    def __lt__(self, other):
        return self.var_data.coords < other.var_data.coords

import logging
from .seqlib import SeqLib


class IdOnlySeqLib(SeqLib):
    """
    Class for counting data with non-variant identifiers and no associated
    FASTQ_ data.
    """

    treeview_class_name = "ID-only SeqLib"

    def __init__(self):
        SeqLib.__init__(self)
        self.identifier_min_count = None
        self.add_label("identifiers")
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from
        a ``.json`` file.
        """
        SeqLib.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )
        try:
            if "min count" in cfg["identifiers"]:
                self.identifier_min_count = int(cfg["identifiers"]["min count"])
            else:
                self.identifier_min_count = 0
        except KeyError as key:
            raise KeyError(
                "Missing required config value {key} [{name}]"
                "".format(key=key, name=self.name)
            )

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for
        dumping to a config file.
        """
        cfg = SeqLib.serialize(self)

        cfg["identifiers"] = dict()
        if self.identifier_min_count > 0:
            cfg["identifiers"]["min count"] = self.identifier_min_count

        return cfg

    def calculate(self):
        """
        Get the identifier counts from the counts file.
        """
        if not self.check_store("/main/identifiers/counts"):
            if self.counts_file is not None:
                self.counts_from_file(self.counts_file)
            else:
                raise ValueError("Missing counts file [{}]".format(self.name))
            self.save_filtered_counts(
                "identifiers", "count >= self.identifier_min_count"
            )

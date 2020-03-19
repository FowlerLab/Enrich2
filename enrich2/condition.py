from .storemanager import StoreManager
from .selection import Selection


class Condition(StoreManager):
    """
    Dummy class for experimental conditions within an 
    :py:class:`~experiment.Experiment`. Required for proper GUI behavior.
    """

    has_store = False  # don't create an HDF5 for Conditions
    treeview_class_name = "Condition"

    def __init__(self):
        StoreManager.__init__(self)
        self.selections = list()

    def configure(self, cfg, configure_children=True):
        StoreManager.configure(self, cfg)
        if configure_children:
            if "selections" not in cfg:
                raise KeyError(
                    "Missing required config value {} [{}]".format(
                        "selections", self.name
                    )
                )

            for sel_cfg in cfg["selections"]:
                sel = Selection()
                sel.configure(sel_cfg)
                self.add_child(sel)

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = StoreManager.serialize(self)
        cfg["selections"] = [child.serialize() for child in self.children]
        return cfg

    def validate(self):
        """
        Calls validate on all child Selections.
        """
        for child in self.children:
            child.validate()

    def _children(self):
        """
        Method bound to the ``children`` property. Returns a list of all 
        :py:class:`~selection.Selection` objects belonging to this object, 
        sorted by name.
        """
        return sorted(self.selections, key=lambda x: x.name)

    def add_child(self, child):
        """
        Add a :py:class:`~selection.Selection`.
        """
        if child.name in self.child_names():
            raise ValueError(
                "Non-unique selection name '{}' [{}]".format(child.name, self.name)
            )
        child.parent = self
        self.selections.append(child)

    def remove_child_id(self, tree_id):
        """
        Remove the reference to a :py:class:`~selection.Selection` with 
        Treeview id *tree_id*.
        """
        self.selections = [x for x in self.selections if x.treeview_id != tree_id]

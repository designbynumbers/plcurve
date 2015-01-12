import spherogram
from plink import LinkEditor
from .. import pdcode as pdc
import Tkinter

def pd_edit(pd=None):
    """pd_edit(pd) -> new PlanarDiagram

    Open a plink editor for the input planar diagram object. On
    completion, returns a new PlanarDiagram object.
    """
    if pd is not None:
        editor = pd.as_spherogram.view()#spherogram.Link(pd.pdcode()).view()
    else:
        editor = LinkEditor()
    Tkinter.mainloop()
    return pdc.PlanarDiagram.from_plink(editor)

if __name__ == "__main__":
    tref = pdc.PlanarDiagram.torus_knot(2,3)
    print repr(pd_edit(tref))

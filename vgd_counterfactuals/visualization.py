import typing as t

import matplotlib.pyplot as plt
import visual_graph_datasets.typing as tv
from matplotlib.backends.backend_pdf import PdfPages
from imageio.v2 import imread


def create_counterfactual_pdf(counterfactual_elements: t.List[tv.VgdElementDict],
                              output_path: str,
                              original_element: t.Optional[tv.VgdElementDict] = None,
                              original_label: t.Optional[str] = None,
                              counterfactual_labels: t.Optional[t.List[str]] = None,
                              counterfactual_color: t.Any = '#FFDCDC',
                              base_fig_size: int = 10,
                              ) -> None:
    """
    Creates a PDF document that visualizes the given ``counterfactual_elements``, where each element is
    displayed on a separate page.

    :param counterfactual_elements: A list of VGD element dicts that represent the counterfactuals to be
        visualized. It is important that these VGD element dicts contain the "image_path" key referencing
        an existing path to a visualization image.
    :param output_path: The path at which to save the PDF file
    :param original_element: Optionally a VGD element dict for the original element from which the
        counterfactuals were created. If this is given an additional page will be added to the beginning
        displaying the original element as well.
    :param original_label: Optionally an additional string, which will be added to the title of the
        original element's page.
    :param counterfactual_labels: Optionally a list of strings, where each string will be used to add to
        the end of that element's title in the PDF.
    :param counterfactual_color: A valid matplotlib color, which will be used the background color for all
        the counterfactual pages to visually separate them from the original element.
    :param base_fig_size: The integer matplotlib fig_size to be used

    :returns: None
    """
    fig_size = (base_fig_size, base_fig_size)

    with PdfPages(output_path) as pdf:

        # The first page we want to be the original prediction of this
        if original_element:
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=fig_size)

            image = imread(original_element['image_path'])
            ax.imshow(image)

            title = (f'ORIGINAL ELEMENT\n'
                     f'{original_element["metadata"]["name"]}')
            if original_label:
                title += f'\n{original_label}'

            fig.suptitle(title)

            pdf.savefig(fig)
            plt.close(fig)

        # Then we need to do all the counterfactuals.
        for index, element in enumerate(counterfactual_elements):
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=fig_size)
            fig.set_facecolor(counterfactual_color)

            image = imread(element['image_path'])
            ax.imshow(image)

            title = (f'COUNTERFACTUAL #{index+1}\n'
                     f'{element["metadata"]["name"]}\n'
                     f'Distance: {element["metadata"]["distance"]:.3f}')
            if counterfactual_labels:
                label = counterfactual_labels[index]
                title += f'\n{label}'

            fig.suptitle(title)

            pdf.savefig(fig)
            plt.close(fig)

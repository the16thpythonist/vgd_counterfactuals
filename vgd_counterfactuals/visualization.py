import typing as t

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from imageio.v2 import imread


def create_counterfactual_pdf(counterfactual_elements: t.List[dict],
                              output_path: str,
                              original_element: t.Optional[dict] = None,
                              original_label: t.Optional[str] = None,
                              counterfactual_labels: t.Optional[t.List[str]] = None,
                              counterfactual_color: t.Any = '#FFDCDC',
                              base_fig_size: int = 10,
                              ) -> None:
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

#!/bin/bash -x

# write to file
bootstrap_template()
{
    pandoc --standalone --toc --template=elegant_bootstrap_menu.html tutorial.md -o index.html --metadata title="Population: new generation population dynamics model with a dynamic population structure"
}

# execute it
bootstrap_template
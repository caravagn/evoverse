install.packages("badger")

badger::badge_custom("Email",
                     "gcaravagn@gmail.com",
                     "blue",
                     "mailto:gcaravagn@gmail.com")

badger::badge_custom("Twitter",
                     "@gcaravagna",
                     "orange",
                     "https://twitter.com/gcaravagna")

badger::badge_custom("Personal webpage",
                     "https://bit.ly/2kc9E6Y",
                     "orangered3",
                     "https://sites.google.com/site/giuliocaravagna/")


badger::badge_custom("Available",
                     "Pipelines",
                     "yellow",
                     "https://caravagn.github.io/evoverse/articles/pipelines.html") %>% cat

badger::badge_custom("Available",
                     "Packages",
                     "yellow",
                     "https://caravagn.github.io/evoverse/articles/packages.html") %>% cat

badger::badge_custom("GitHub Pages",
                     "https://caravagn.github.io/evoverse/",
                     "yellow",
                     "https://caravagn.github.io/evoverse") %>% cat

badger::badge_custom("Part of",
                     "evoverse",
                     "blue",
                     "https://caravagn.github.io/evoverse") %>% cat

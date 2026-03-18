This branch is only here to publish documentation on github.io server.

Follow these steps to do it:

- Copy the generated html directory of each version in `docs/` folder
- Create a file named `.nojekyll` in each version folder. Otherwise, javascript rendering will be deactivated
- Create or update the `docs/index.html` file to make it redirect to the latest doc version, for example :
  ```html
  <head>
    <meta http-equiv='refresh' content='0; URL=1.2.0/index.html'>
  </head>
  ```
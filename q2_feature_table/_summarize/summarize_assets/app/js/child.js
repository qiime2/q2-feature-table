document.addEventListener('DOMContentLoaded', function() {
  parent.postMessage(document.body.scrollHeight, '*');
});

<script>
   function showPage(pageId) {
      document.querySelectorAll(".page").forEach(page => {
         page.classList.remove('active');
      });
      
      document.getElementById(pageId).classList.add('active');
      
      window.scrollTo({
         top: 0,
         behavior: 'instant'
      });
   }
</script>
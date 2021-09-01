import json
from django.contrib import messages
from django.urls import reverse
from django.views.generic import FormView, TemplateView
from cluster import forms


class ConfigCreate(FormView):
    form_class = forms.ConfigForm
    template_name = "cluster/config_create.html"

    def form_valid(self, form):
        # Store the values in the form as a message as Json.
        messages.add_message(self.request,
                             messages.INFO,
                             json.dumps(form.cleaned_data))

        return super().form_valid(form)

    def get_success_url(self):
        return reverse("cluster:detail")


class ConfigDetail(TemplateView):
    template_name = "cluster/config_detail.html"
    content_type = "text/plain; charset=utf-8"

    def get_context_data(self, **kwargs):
        kwargs = super().get_context_data(**kwargs)

        # Retrieve Json data from the messages.
        for message in messages.get_messages(self.request):
            kwargs.update(json.loads(message.message))

        return kwargs

    def render_to_response(self, context, **response_kwargs):
        response = super().render_to_response(context, **response_kwargs)

        # Set the reponse type as a downloadable file, instead a browser view.
        response.headers.setdefault(
            "Content-Disposition", "attachment; filename=config.py")

        return response

from django.core.exceptions import ValidationError


def validate_cysbond(values: str) -> None:
    """Validate that value can be translated to a list of CYS bonds."""
    pairs = values.split(",")
    for cysBond in pairs:
        try:
            (cysA, cysB) = cysBond.split(":")
        except ValueError:
            raise ValidationError(
                f"Invalid value detected ({cysBond}): must be integer:integer")
        try:
            if int(cysA) < 1 or int(cysB) < 1:
                raise ValidationError(
                    f"Invalid value detected ({cysA}:{cysB}): less than 1")
        except ValueError:
            raise ValidationError(
                f"Invalid value detected ({cysA}:{cysB}): not an integer")
        if int(cysA) > int(cysB):
            raise ValidationError(
                f"Invalid value detected ({cysA}:{cysB}): {cysB} should be "
                f"higher than {cysA}")

def validate_temps(values: str) -> None:
    """Validate the temperatures."""
    temps = values.split(",")
    for temp in temps:
        try:
            tempInt = int(temp)
        except ValueError:
            raise ValidationError(
                f"Invalid value detected ({temp}): must be an integer")
        if tempInt < 0:
            raise ValidationError(
                f"Invalid value detected ({temp}): ºK temperatures are always above 0")


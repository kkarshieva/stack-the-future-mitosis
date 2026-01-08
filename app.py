from flask import Flask, render_template

app = Flask(__name__)


@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")


@app.route("/onboarding", methods=["GET"])
def onboarding():
    return render_template("onboarding.html")


@app.route("/matching", methods=["GET"])
def matching():
    return render_template("matching.html")


@app.route("/dashboard", methods=["GET"])
def dashboard():
    return render_template("dashboard.html")
